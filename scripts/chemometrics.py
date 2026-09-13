"""
chemometrics.py  -  multivariate calibration (PLS / PCR / PCA) for spectra, honestly.

This is the FOURTH family (after spectra, calibration, crystal). Univariate
`calibration.py` regresses ONE band area on concentration; this regresses the WHOLE
spectrum, which is what you reach for when no single band is selective (overlap,
matrix effects, scatter/loading variation). The figures are the standard PLS read-out:
RMSECV-vs-components, the regression-coefficient spectrum, scores, and loadings.

The correctness traps this guards against (the reason it lives in THIS skill and not
a notebook):

  * CV LEAKAGE. Preprocessing that uses cross-sample statistics - MSC's reference
    spectrum, mean-centring, autoscaling - must be fit on the TRAINING rows of each
    fold only. Fit it once on the whole set and the test rows have leaked into
    training, so RMSECV is optimistic. `cross_validate` re-fits a FRESH `Preprocessor`
    inside every fold. (SNV and derivatives are per-sample, so they don't leak - but
    the machinery is uniform so you can't get it wrong by accident.)

  * OVER-FITTING THE COMPONENT COUNT. The global-minimum RMSECV almost always sits at
    too many latent variables. `choose_n_components` applies a PARSIMONY rule: the
    fewest components whose RMSECV is within `cfg.lv_parsimony_tol` of the global min.

  * THE OPTIMISTIC-CV TRAP. Leave-one-out over replicate spectra that share a
    concentration LEVEL answers "predict a new sample at a level I've already seen",
    not "predict a new level". The two give different errors. Every CV result NAMES
    which question its scheme answers, and `cross_validate` WARNs if you LOO over data
    whose y replicates (telling you to pass `groups=` for leave-one-level-out).

scikit-learn is a LAZY dependency: it is imported only when a model is actually fit,
mirroring the crystal family's heavy deps. An FTIR or univariate-calibration job never
imports it.

    Preprocessor(steps, cfg)              leakage-aware SNV / MSC / d1 / d2 / centre / autoscale
    pls_fit(X, y, nc, cfg, pre=...)       -> model dict (coef, scores, loadings, weights, yhat...)
    pcr_fit / pca_fit                     PCR regression / PCA decomposition (same dict shape)
    cross_validate(X, y, nc, cfg, ...)    leakage-safe RMSECV + held-out predictions, scheme NAMED
    component_scan(X, y, cfg, ...)        RMSECV vs components, one curve per preprocessing
                                          (one fit per fold at the cap, truncated to 1..cap)
    choose_n_components(scan, cfg)        the parsimony pick
    vip(model)                            VIP scores (which wavenumbers drive the model)
    plot_rmsecv / plot_coefficients / plot_scores / plot_loadings   ax-level panels
    diagnostics_figure(X, y, x, cfg, ...) the gated + QA'd 2x2 composite deliverable
    methods_text(model, cv, choice, cfg)  paste-ready methods paragraph
    permutation_test / corrected_paired_t / rpd / rpiq   the validation trio
                                          (per-fold preprocessing cached across permutations)

(methods_text, not methods_report — the latter is report.py's name, and the bundler
flattens every module into one namespace, so module-level names must be unique.)
"""
from __future__ import annotations
import functools
import numpy as np
from . import style, verify

_SKLEARN_HINT = ("Required package not found: scikit-learn. Install with:\n"
                 "    python -m pip install scikit-learn\n"
                 "(the chemometrics family needs sklearn for PLS / PCA / PCR).")


def _pls_cls():
    try:
        from sklearn.cross_decomposition import PLSRegression
        return PLSRegression
    except ImportError as e:
        raise ImportError(_SKLEARN_HINT) from e


def _pca_cls():
    """PCA with the deterministic full SVD pinned ONCE, so the CV scan, the deployed pca_fit /
    pcr_fit model and the nested truncation all decompose the same way (sklearn's 'auto'
    policy may pick a randomised solver on wide data)."""
    try:
        from sklearn.decomposition import PCA
        return functools.partial(PCA, svd_solver="full")
    except ImportError as e:
        raise ImportError(_SKLEARN_HINT) from e


def _linreg_cls():
    try:
        from sklearn.linear_model import LinearRegression
        return LinearRegression
    except ImportError as e:
        raise ImportError(_SKLEARN_HINT) from e


# ===================================================== leakage-aware preprocessing
_STATELESS = {"none", "snv", "d1", "d2"}
_STATEFUL = {"msc", "center", "autoscale"}


class Preprocessor:
    """A small preprocessing pipeline whose STATEFUL steps are fit on training rows
    only, so it is safe to drive cross-validation with (the whole reason it exists).

    steps: a list like ["snv"] or ["d1", "center"], or a "+"-joined string
    ("snv+center"). Order is applied left to right.
      stateless (identical train/test, no leakage):
        "none"  pass-through
        "snv"   standard normal variate - per-row centre & scale (scatter/path-length)
        "d1"    Savitzky-Golay 1st derivative (window from cfg.deriv_smooth)
        "d2"    Savitzky-Golay 2nd derivative (resolves overlap; amplifies noise)
      stateful (parameters estimated from TRAINING rows -> must go through fit):
        "msc"       multiplicative scatter correction against the training-mean spectrum
        "center"    subtract the training per-wavenumber mean
        "autoscale" centre and divide by the training per-wavenumber SD

    Use fit_transform on a training block and transform on the matching test block;
    cross_validate / component_scan do this for every fold automatically.

    Note on the scale convention: SNV and `autoscale` use the POPULATION standard
    deviation (ddof=0), matching sklearn's StandardScaler and the `chemotools` library.
    R/prospectr use the sample SD (ddof=1) -- a constant sqrt(n/(n-1)) factor. The ddof=0
    choice is intentional and is locked by the chemometrics validation suite
    (skill_validation/chemometrics/)."""

    def __init__(self, steps, cfg=None):
        if isinstance(steps, str):
            steps = [s.strip() for s in steps.split("+") if s.strip()]
        steps = list(steps) if steps is not None and len(steps) else ["none"]
        bad = [s for s in steps if s not in (_STATELESS | _STATEFUL)]
        if bad:
            raise ValueError(f"unknown preprocessing step(s): {bad} "
                             f"(allowed: {sorted(_STATELESS | _STATEFUL)})")
        self.steps = steps
        self.cfg = cfg
        self._state = {}        # per-step-index state for the stateful steps

    @property
    def name(self):
        return "+".join(self.steps)

    def _deriv(self, X, order):
        from scipy.signal import savgol_filter
        win, poly = (self.cfg.deriv_smooth if self.cfg is not None else (17, 3))
        w = win if win % 2 == 1 else win + 1
        w = min(w, X.shape[1] - (1 - X.shape[1] % 2))   # window can't exceed n_features
        if w <= poly:
            import warnings
            warnings.warn(
                f"Savitzky-Golay d{order} skipped: usable window ({w}) <= polyorder "
                f"({poly}) for {X.shape[1]} features; returning the input UNCHANGED "
                f"(no derivative applied).", stacklevel=2)
            return X
        return savgol_filter(X, w, poly, deriv=order, axis=1)

    def _run(self, X, fit):
        X = np.asarray(X, float)
        for i, s in enumerate(self.steps):
            st = self._state.get(i, {})
            if s == "none":
                pass
            elif s == "snv":
                mu = X.mean(axis=1, keepdims=True)
                sd = X.std(axis=1, keepdims=True)
                X = (X - mu) / np.where(sd == 0, 1.0, sd)
            elif s == "d1":
                X = self._deriv(X, 1)
            elif s == "d2":
                X = self._deriv(X, 2)
            elif s == "msc":
                if fit:
                    st["ref"] = X.mean(axis=0)
                ref = st["ref"]
                out = np.empty_like(X)
                for r in range(X.shape[0]):
                    b1, b0 = np.polyfit(ref, X[r], 1)        # X[r] ~ b1*ref + b0
                    out[r] = (X[r] - b0) / (b1 if b1 != 0 else 1.0)
                X = out
            elif s == "center":
                if fit:
                    st["mean"] = X.mean(axis=0)
                X = X - st["mean"]
            elif s == "autoscale":
                if fit:
                    st["mean"] = X.mean(axis=0)
                    sd = X.std(axis=0)
                    st["std"] = np.where(sd == 0, 1.0, sd)
                X = (X - st["mean"]) / st["std"]
            if fit:
                self._state[i] = st
        return X

    def fit(self, X):
        self._run(X, fit=True)
        return self

    def transform(self, X):
        return self._run(X, fit=False)

    def fit_transform(self, X):
        return self._run(X, fit=True)


# ===================================================== EMSC
def emsc(X, wavenumbers, reference=None, poly_order=2, interferents=None,
         weights=None, return_model=False):
    """Extended Multiplicative Signal Correction (Martens & Stark 1991; Afseth & Kohler 2012).

    Fits each spectrum    x ~= b*ref + sum_k a_k * z^k  (+ sum_j d_j * interferent_j),
    where z = wavenumbers affine-scaled to [-1, 1], then returns the corrected spectrum
        x_corr = (x - baseline - interferents) / b.
    `reference` defaults to the MEAN spectrum (the EMSC convention). Validated bit-for-bit
    (~5e-16) against the biospectools / Kohler MATLAB gold standard.

    The principled upgrade over SNV/MSC when the problem is multiplicative loading plus a
    baseline. `interferents` (rows of an array, e.g. a pure-excipient spectrum) models known
    constituents and subtracts them — the project's "model the matrix away" move.

    Caveat on b: with reference=mean AND a baseline term, b is NOT a physical dilution factor
    (the constant polynomial column absorbs part of the offset). To recover a known multiplicative
    factor, use a PURE-component reference with poly_order=0.

    Leakage note: a FIXED (pure/external) reference + interferents is leak-free — apply it before
    CV. The mean-reference variant is cross-sample; refit it per fold if used inside cross_validate.

    Returns the corrected array (same leading dim as X); also the coef/model dict if return_model.
    """
    Xin = np.asarray(X, float)
    X2 = np.atleast_2d(Xin)
    wn = np.asarray(wavenumbers, float)
    ref = X2.mean(0) if reference is None else np.asarray(reference, float)
    rng = wn.max() - wn.min()
    z = 2.0 * (wn - wn.min()) / (rng if rng else 1.0) - 1.0
    polys = np.vstack([z ** k for k in range(poly_order + 1)]).T        # (nwn, p+1)
    blocks = [ref[:, None], polys]
    I = None
    if interferents is not None:
        I = np.atleast_2d(np.asarray(interferents, float))             # (n_interf, nwn)
        blocks.append(I.T)
    M = np.hstack(blocks)                                              # [ ref | polys | interferents ]
    npoly = poly_order + 1
    w = None if weights is None else np.sqrt(np.asarray(weights, float))
    Mfit = M if w is None else M * w[:, None]

    corrected = np.empty_like(X2)
    coefs = np.empty((X2.shape[0], M.shape[1]))
    for i in range(X2.shape[0]):
        target = X2[i] if w is None else X2[i] * w
        c, *_ = np.linalg.lstsq(Mfit, target, rcond=None)
        b = c[0]
        baseline = polys @ c[1:1 + npoly]
        interf = (I.T @ c[1 + npoly:]) if I is not None else 0.0
        corrected[i] = (X2[i] - baseline - interf) / (b if b != 0 else 1.0)
        coefs[i] = c
    out = corrected if Xin.ndim == 2 else corrected[0]
    if return_model:
        return out, {"coefs": coefs, "ref": ref, "M": M, "poly_order": poly_order,
                     "n_interferents": 0 if I is None else I.shape[0],
                     "residual": X2 - coefs @ M.T}
    return out


# ===================================================== input gate
def _check_xy(X, y, cfg, n_components=None):
    X = np.asarray(X, float)
    y = np.asarray(y, float).ravel()
    f = []
    if X.ndim != 2:
        f.append(("FAIL", f"X must be 2D (n_samples, n_features); got shape {X.shape}"))
        return verify._resolve(f, cfg, "chemometrics"), X, y
    if X.shape[0] != y.size:
        f.append(("FAIL", f"X has {X.shape[0]} rows but y has {y.size}"))
    if not np.isfinite(X).all():
        f.append(("FAIL", "X has non-finite values (NaN/inf) - check the ingest/baseline"))
    if not np.isfinite(y).all():
        f.append(("FAIL", "y has non-finite values"))
    if X.shape[0] < 4:
        f.append(("WARN", f"only {X.shape[0]} samples - CV and component estimates are unstable"))
    if n_components is not None and n_components > X.shape[0] - 1:
        f.append(("WARN", f"n_components={n_components} >= n_samples ({X.shape[0]}); over-fit risk"))
    if not f:
        f.append(("INFO", f"X {X.shape}, y n={y.size} - OK"))
    verify._resolve(f, cfg, "chemometrics")
    return f, X, y


# ===================================================== model fits (on the FULL set)
def _coef_vector(model, n_features):
    coef = np.ravel(np.asarray(model.coef_))
    if coef.size != n_features:                 # sklearn version differences in orientation
        coef = np.ravel(np.asarray(model.coef_).T)
    return coef


def _fit_stats(y, yhat):
    resid = y - yhat
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - float(np.sum(resid ** 2)) / ss_tot if ss_tot else float("nan")
    rmse = float(np.sqrt(np.mean(resid ** 2)))
    return resid, float(r2), rmse


def pls_fit(X, y, n_components, cfg, pre=None):
    """Fit PLS regression on the WHOLE set. Use this for the interpretable read-out
    (coefficients, scores, loadings) and the calibration fit; use cross_validate for
    the HONEST error. Returns a model dict."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, y = _check_xy(X, y, cfg, n_components)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    nc = min(n_components, X.shape[0] - 1, Xp.shape[1])
    PLS = _pls_cls()
    m = PLS(n_components=nc, scale=cfg.pls_scale).fit(Xp, y)
    coef = _coef_vector(m, Xp.shape[1])
    yhat = np.ravel(m.predict(Xp))
    resid, r2, rmsec = _fit_stats(y, yhat)
    return {"kind": "pls", "model": m, "pre": pp, "pre_name": pp.name,
            "n_components": int(nc), "coef": coef,
            "scores": np.asarray(m.x_scores_), "loadings": np.asarray(m.x_loadings_),
            "weights": np.asarray(m.x_weights_), "y_loadings": np.ravel(m.y_loadings_),
            "yhat": yhat, "resid": resid, "r2_cal": r2, "rmsec": rmsec, "y": y,
            "caption": f"PLS ({pp.name}, {nc} LV): R2(cal)={r2:.3f}, RMSEC={rmsec:.3g}"}


def pca_fit(X, cfg, n_components=None, pre=None):
    """PCA decomposition for exploratory scores/loadings (no y). sklearn PCA centres
    internally, so `pre` need not include 'center'."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, _y = _check_xy(X, np.zeros(X.shape[0]), cfg)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    nc = n_components or min(cfg.pls_max_components, *Xp.shape)
    nc = int(min(nc, *Xp.shape))            # cap at min(n_samples, n_features) (sklearn's own limit)
    PCA = _pca_cls()
    m = PCA(n_components=nc).fit(Xp)
    return {"kind": "pca", "model": m, "pre": pp, "pre_name": pp.name,
            "n_components": int(nc), "scores": np.asarray(m.transform(Xp)),
            "loadings": np.asarray(m.components_.T),
            "explained": np.asarray(m.explained_variance_ratio_),
            "caption": f"PCA ({pp.name}, {nc} PC): "
                       f"{100*float(np.sum(m.explained_variance_ratio_)):.1f}% variance"}


def pcr_fit(X, y, n_components, cfg, pre=None):
    """Principal-component regression: PCA, then OLS of y on the first n_components
    scores. Coefficients are mapped back to wavenumber space for the spectrum plot."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, y = _check_xy(X, y, cfg, n_components)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    nc = min(n_components, X.shape[0] - 1, Xp.shape[1])
    PCA = _pca_cls()
    pca = PCA(n_components=nc).fit(Xp)
    T = pca.transform(Xp)
    reg = _linreg_cls()().fit(T, y)
    yhat = np.ravel(reg.predict(T))
    coef = pca.components_.T @ np.ravel(reg.coef_)     # back to X space (p,)
    resid, r2, rmsec = _fit_stats(y, yhat)
    return {"kind": "pcr", "model": (pca, reg), "pre": pp, "pre_name": pp.name,
            "n_components": int(nc), "coef": coef, "scores": np.asarray(T),
            "loadings": np.asarray(pca.components_.T),
            "explained": np.asarray(pca.explained_variance_ratio_),
            "yhat": yhat, "resid": resid, "r2_cal": r2, "rmsec": rmsec, "y": y,
            "caption": f"PCR ({pp.name}, {nc} PC): R2(cal)={r2:.3f}, RMSEC={rmsec:.3g}"}


# ===================================================== cross-validation (leakage-safe)
def _make_folds(n, y, groups, scheme, cfg):
    """Return (folds, scheme_name, question). folds = list of (train_idx, test_idx)."""
    idx = np.arange(n)
    if groups is not None:
        groups = np.asarray(groups)
        uniq = list(dict.fromkeys(groups.tolist()))     # stable order
        folds = [(idx[groups != g], idx[groups == g]) for g in uniq]
        return (folds, f"leave-one-group-out ({len(uniq)} groups)",
                "predict a NEW group (e.g. an unseen concentration level)")
    sch = _scheme(scheme, cfg)
    if sch in ("auto", "loo"):
        folds = [(idx[idx != i], np.array([i])) for i in idx]
        return (folds, "leave-one-out (per sample)",
                "predict a new sample at a SEEN level "
                "(optimistic for a new level when samples replicate levels)")
    if sch == "kfold":
        # stratified by LEVEL: the replicates of one y value stay together (holding them out
        # one at a time leaves their level in training - leave-one-out's optimism in disguise),
        # and the sorted levels are interleaved across folds so every fold spans the response
        # range (contiguous blocks would make the end folds pure extrapolation).
        levels = np.unique(y)
        k = int(min(cfg.cv_folds, len(levels)))
        parts = [idx[np.isin(y, levels[i::k])] for i in range(k)]
        parts = [p for p in parts if len(p)]
        folds = [(np.setdiff1d(idx, p), p) for p in parts]
        return (folds, f"{k}-fold (stratified by level; a level's replicates stay together)",
                "predict held-out LEVELS (each fold spans the y range)")
    raise ValueError(f"unknown cv_scheme: {sch!r} (use 'auto'/'loo'/'kfold' or pass groups=)")


def _fit_predict(Xtr, ytr, Xte, nc, cfg, method):
    """One model at `nc` components, from scratch (the reference path; `_fit_predict_nested`
    is the fast one the scan uses)."""
    nc = int(min(nc, len(Xtr) - 1, Xtr.shape[1]))
    if method == "pls":
        m = _pls_cls()(n_components=nc, scale=cfg.pls_scale).fit(Xtr, ytr)
        return np.ravel(m.predict(Xte))
    if method == "pcr":
        pca = _pca_cls()(n_components=nc).fit(Xtr)
        reg = _linreg_cls()().fit(pca.transform(Xtr), ytr)
        return np.ravel(reg.predict(pca.transform(Xte)))
    raise ValueError(f"unknown method: {method!r} (use 'pls' or 'pcr')")


def _fit_predict_nested(Xtr, ytr, Xte, cap, cfg, method):
    """ONE fit at `cap` components, predictions for a = 1..cap by truncation: returns an
    array (len(Xte), cap) whose column a-1 equals `_fit_predict(..., a, ...)`.

    PLS1 (NIPALS): the first a weights/loadings of a cap-component model ARE the a-component
    model's, and P'W is upper triangular, so the leading block of the rotation matrix
    R = W (P'W)^-1 is the a-component rotation. Hence the a-component prediction is
    intercept + s * T[:, :a] @ q[:a] with T = model.transform(X), q the y-loadings and s the
    y scale sklearn applied (its ddof=1 SD when scale=True, 1 otherwise). s is read off the
    fitted model itself - the least-squares ratio of (predict - intercept) to T @ q on the
    training rows - so no sklearn convention is hand-replicated. PCR: PCA scores are
    orthogonal and mean-zero, so the OLS coefficients on the first a scores are the leading a
    coefficients of the cap-score regression. Both identities are pinned column by column in
    tests/test_nested_scan.py; the cap column is also checked here against the model's own
    predict() and a mismatch RAISES (a silent fallback would let a convention change in
    scikit-learn switch the algorithm behind a print)."""
    cap = int(min(cap, len(Xtr) - 1, Xtr.shape[1]))
    if method == "pls":
        m = _pls_cls()(n_components=cap, scale=cfg.pls_scale).fit(Xtr, ytr)
        q = np.ravel(m.y_loadings_)                            # (cap,)
        y_mean = float(np.ravel(m.intercept_)[0])
        t_tr = np.asarray(m.transform(Xtr)) @ q
        p_tr = np.ravel(m.predict(Xtr)) - y_mean
        tt = float(t_tr @ t_tr)
        scale = float(p_tr @ t_tr) / tt if tt > 0 else 1.0
        out = y_mean + scale * np.cumsum(np.asarray(m.transform(Xte)) * q, axis=1)
        check = np.ravel(m.predict(Xte))
    elif method == "pcr":
        pca = _pca_cls()(n_components=cap).fit(Xtr)
        Ttr, Tte = pca.transform(Xtr), pca.transform(Xte)
        reg = _linreg_cls()().fit(Ttr, ytr)
        b = np.ravel(reg.coef_)
        out = float(np.ravel(reg.intercept_)[0]) + np.cumsum(Tte * b, axis=1)
        check = np.ravel(reg.predict(Tte))
    else:
        raise ValueError(f"unknown method: {method!r} (use 'pls' or 'pcr')")
    tol = 1e-8 * max(float(np.max(np.abs(check))), 1.0)
    if not np.allclose(out[:, -1], check, rtol=0, atol=tol):
        raise RuntimeError(
            f"chemometrics: nested {method} truncation disagrees with predict() at {cap} "
            f"components (max |d| = {np.max(np.abs(out[:, -1] - check)):.2e}); the truncation "
            f"identity no longer holds for this scikit-learn - use _fit_predict per count")
    return out


def _fold_blocks(X, cfg, pre, folds):
    """Preprocess each fold ONCE: [(tr, te, Xtr_p, Xte_p)], a FRESH Preprocessor fit on the
    training rows only. Every step is a function of X alone (SNV/derivatives per row; MSC,
    centre, autoscale from the training rows), so the blocks are reusable across component
    counts and across y-permutations. A pipeline of row-wise (stateless) steps only is applied
    to X once and sliced - each row's transform does not depend on which rows share the fold."""
    pp = Preprocessor(pre, cfg)
    if all(step in _STATELESS for step in pp.steps):
        Xp = pp.fit_transform(X)
        return [(tr, te, Xp[tr], Xp[te]) for tr, te in folds]
    blocks = []
    for tr, te in folds:
        pp = Preprocessor(pre, cfg)                      # FRESH per fold
        Xtr = pp.fit_transform(X[tr])                    # fit on TRAIN rows only
        Xte = pp.transform(X[te])
        blocks.append((tr, te, Xtr, Xte))
    return blocks


def _cv_predict(X, y, nc, cfg, pre, folds, method, blocks=None):
    """Held-out predictions, preprocessing re-fit inside every fold (no leakage)."""
    if blocks is None:
        blocks = _fold_blocks(X, cfg, pre, folds)
    pred = np.full(len(y), np.nan)
    for tr, te, Xtr, Xte in blocks:
        pred[te] = _fit_predict(Xtr, y[tr], Xte, nc, cfg, method)
    return pred


def _cv_predict_scan(X, y, cap, cfg, pre, folds, method):
    """Held-out predictions for EVERY component count 1..cap at once: array (n, cap), one
    fit per fold at `cap`, truncated (see _fit_predict_nested)."""
    pred = np.full((len(y), cap), np.nan)
    for tr, te, Xtr, Xte in _fold_blocks(X, cfg, pre, folds):
        p = _fit_predict_nested(Xtr, y[tr], Xte, cap, cfg, method)
        assert p.shape[1] == cap          # component_scan caps at min_train-1 and n_features
        pred[te] = p
    return pred


def _scheme(scheme, cfg):
    return (scheme or cfg.cv_scheme).lower()


def _cv_stats(y, pred):
    """(resid, rmsecv, r2cv, bias) of held-out predictions; NaN-aware."""
    resid = y - pred
    rmsecv = float(np.sqrt(np.nanmean(resid ** 2)))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2cv = 1 - float(np.nansum(resid ** 2)) / ss_tot if ss_tot else float("nan")
    bias = float(np.nanmean(pred - y))
    return resid, rmsecv, r2cv, bias


def _optimistic_loo(groups, scheme, y, cfg):
    return (groups is None and _scheme(scheme, cfg) in ("auto", "loo")
            and np.unique(y).size < y.size)


def cross_validate(X, y, n_components, cfg, pre=None, groups=None, scheme=None,
                   method="pls", report=True):
    """Leakage-safe cross-validation. Re-fits the preprocessor inside each fold and
    NAMES which prediction question the scheme answers. Pass groups=level-labels for
    leave-one-level-out. Returns a dict with pred, rmsecv, r2cv, bias, scheme, question."""
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, y = _check_xy(X, y, cfg, n_components)
    folds, scheme_name, question = _make_folds(len(y), y, groups, scheme, cfg)
    pred = _cv_predict(X, y, n_components, cfg, pre, folds, method)
    resid, rmsecv, r2cv, bias = _cv_stats(y, pred)
    if report:
        f = []
        if _optimistic_loo(groups, scheme, y, cfg):
            f.append(("WARN", "LOO over data with replicated y answers 'predict a SEEN level'; "
                              "for a NEW level pass groups=level-labels (leave-one-level-out)"))
        f.append(("INFO", f"{scheme_name}: RMSECV={rmsecv:.3g} {cfg.conc_unit}, "
                          f"R2cv={r2cv:.3f}, bias={bias:+.3g}, {len(folds)} folds"))
        verify._resolve(f, cfg, "cross-validate")
    return {"pred": pred, "resid": resid, "rmsecv": rmsecv, "r2cv": r2cv, "bias": bias,
            "n_components": int(n_components), "method": method,
            "pre_name": Preprocessor(pre, cfg).name, "scheme": scheme_name,
            "question": question, "n_splits": len(folds),
            "caption": f"RMSECV={rmsecv:.3g} {cfg.conc_unit} "
                       f"({scheme_name}); answers: {question}"}


def component_scan(X, y, cfg, max_components=None, pre_options=None, groups=None,
                   scheme=None, method="pls"):
    """RMSECV against component count, ONE curve per preprocessing, all on the same
    folds (so the curves are comparable). Drives choose_n_components and plot_rmsecv.
    One model per fold at the cap, truncated to every smaller count (_fit_predict_nested)."""
    _, X, y = _check_xy(X, y, cfg)
    n = len(y)
    folds, scheme_name, question = _make_folds(n, y, groups, scheme, cfg)
    min_train = min(len(tr) for tr, _ in folds)
    cap = int(min(max_components or cfg.pls_max_components,
                  min_train - 1, X.shape[1]))
    cap = max(cap, 1)
    if pre_options is None:
        pre_options = (cfg.pls_preprocess,)
    comps = list(range(1, cap + 1))
    curves = {}
    f = []
    if _optimistic_loo(groups, scheme, y, cfg):
        f.append(("WARN", "LOO over replicated y is optimistic for a NEW level; "
                          "pass groups=level-labels for leave-one-level-out"))
    for pre in pre_options:
        name = Preprocessor(pre, cfg).name
        pred = _cv_predict_scan(X, y, cap, cfg, pre, folds, method)
        rms = [_cv_stats(y, pred[:, i])[1] for i in range(cap)]
        curves[name] = rms
        f.append(("INFO", f"{name:14s} RMSECV " + " ".join(f"{r:.2f}" for r in rms)))
    verify._resolve(f, cfg, f"component-scan[{scheme_name}]")
    return {"curves": curves, "components": comps, "max": cap, "method": method,
            "scheme": scheme_name, "question": question}


def choose_n_components(scan, cfg, tol=None):
    """Parsimony pick: the FEWEST components whose RMSECV is within `tol` (fraction)
    of the global-min RMSECV across all scanned preprocessings. Ties broken by lower
    RMSECV. Avoids the over-fit you get from taking the global argmin."""
    tol = cfg.lv_parsimony_tol if tol is None else tol
    flat = [(pre, nc, scan["curves"][pre][i])
            for pre in scan["curves"] for i, nc in enumerate(scan["components"])]
    gmin = min(v for _, _, v in flat)
    cands = [t for t in flat if t[2] <= gmin * (1 + tol)]
    pre, nc, rmsecv = min(cands, key=lambda t: (t[1], t[2]))
    return {"pre": pre, "n_components": int(nc), "rmsecv": float(rmsecv),
            "global_min": float(gmin), "tol": tol,
            "caption": f"{pre}, {nc} component(s) (RMSECV={rmsecv:.3g}); parsimony: "
                       f"fewest components within {int(tol*100)}% of the global min {gmin:.3g}"}


def vip(model):
    """VIP scores: per-wavenumber importance for a PLS model (>1 ~ above-average
    contribution). Needs a PLS model dict (weights, scores, y_loadings)."""
    W = np.asarray(model["weights"])                    # (p, A)
    T = np.asarray(model["scores"])                     # (n, A)
    q = np.ravel(model["y_loadings"])                   # (A,)
    p, A = W.shape
    Wn = W / np.where(np.linalg.norm(W, axis=0) == 0, 1.0, np.linalg.norm(W, axis=0))
    ssy = (q ** 2) * np.sum(T ** 2, axis=0)             # explained Y SS per component
    total = float(np.sum(ssy))
    if total == 0:
        return np.zeros(p)
    return np.sqrt(p * ((Wn ** 2) @ ssy) / total)


# ===================================================== ax-level plot helpers
_CYCLE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
_MARKERS = ["o", "s", "^", "D", "v", "P"]


def _xaxis(ax, x, cfg):
    """Honour FTIR inversion / cfg.x_limits for a wavenumber x-axis."""
    if cfg.x_limits:
        ax.set_xlim(*cfg.x_limits)
    elif cfg.domain == "ftir":
        ax.set_xlim(float(np.max(x)), float(np.min(x)))


def plot_rmsecv(ax, scan, cfg, choice=None):
    """RMSECV vs component count, one line per preprocessing; rings the chosen point."""
    comps = scan["components"]
    for i, (name, rms) in enumerate(scan["curves"].items()):
        ax.plot(comps, rms, "-" + _MARKERS[i % len(_MARKERS)],
                color=_CYCLE[i % len(_CYCLE)], ms=4, lw=1, label=name)
    if choice is not None:
        ax.scatter([choice["n_components"]], [choice["rmsecv"]], s=90, zorder=5,
                   facecolor="none", edgecolor="k", linewidth=1.2)
    ax.set_xticks(comps)
    ax.set_xlabel("Components (latent variables)")
    ax.set_ylabel(f"RMSECV / {cfg.conc_unit}")
    if len(scan["curves"]) > 1:
        ax.legend(fontsize=6.5)
    return ax


def plot_coefficients(ax, model, x, cfg, ref=None, ref_label="reference", highlight=None):
    """Regression-coefficient spectrum (which wavenumbers drive the prediction), with
    an optional pure-component reference overlay (scaled) and a highlighted band."""
    x = np.asarray(x, float)
    coef = np.asarray(model["coef"], float)
    ax.axhline(0, color="0.7", lw=0.5)
    if highlight is not None:
        ax.axvspan(min(highlight), max(highlight), color="0.82", alpha=0.6, zorder=0)
    if ref is not None:
        ref = np.asarray(ref, float)
        scale = np.max(np.abs(coef)) / (np.max(np.abs(ref)) or 1.0)
        ax.fill_between(x, ref * scale, color="#56B4E9", alpha=0.35, lw=0, label=ref_label)
    ax.plot(x, coef, color="#222222", lw=0.8, zorder=3, label="coef.")
    _xaxis(ax, x, cfg)
    ax.set_xlabel(r"Wavenumber / cm$^{-1}$" if cfg.domain == "ftir"
                  else r"2$\theta$ / deg")
    ax.set_ylabel("Regression coef. (a.u.)")
    ax.legend(fontsize=6, loc="best")
    return ax


def plot_scores(ax, model, cfg, color_by=None, comps=(0, 1), cbar=None, cbar_label=None):
    """Scores scatter (the calibration's sample structure). Colour by concentration to
    see the gradient. Returns the PathCollection so the caller can add a colorbar
    (do that AFTER finalize on a constrained-layout figure so it doesn't collide)."""
    T = np.asarray(model["scores"])
    if T.ndim == 1:
        T = T[:, None]
    i, j = comps
    one = T.shape[1] < 2 or j is None or j >= T.shape[1]     # a 1-component model: scores vs sample
    xs_ = T[:, i]
    ys_ = np.arange(T.shape[0]) if one else T[:, j]
    kw = dict(s=26, edgecolor="white", linewidth=0.3)
    if color_by is not None:
        sc = ax.scatter(xs_, ys_, c=np.asarray(color_by, float), cmap="viridis", **kw)
    else:
        sc = ax.scatter(xs_, ys_, color=_CYCLE[0], **kw)
    tag = "LV" if model.get("kind") == "pls" else "PC"
    ax.set_xlabel(f"{tag}{i+1} score (a.u.)")
    ax.set_ylabel("sample (index)" if one else f"{tag}{j+1} score (a.u.)")
    if cbar is not None and color_by is not None:
        cb = cbar.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label(cbar_label or cfg.conc_unit, fontsize=7)
    return sc


def plot_loadings(ax, model, x, cfg, comps=(0, 1), highlight=None):
    """Loadings (the latent variables themselves) vs wavenumber."""
    x = np.asarray(x, float)
    P = np.asarray(model["loadings"])
    tag = "LV" if model.get("kind") == "pls" else "PC"
    ax.axhline(0, color="0.7", lw=0.5)
    if highlight is not None:
        ax.axvspan(min(highlight), max(highlight), color="0.82", alpha=0.6, zorder=0)
    for k, c in enumerate(comps):
        ax.plot(x, P[:, c], color=_CYCLE[k % len(_CYCLE)], lw=0.8,
                alpha=1.0 if k == 0 else 0.75, label=f"{tag}{c+1}")
    _xaxis(ax, x, cfg)
    ax.set_xlabel(r"Wavenumber / cm$^{-1}$" if cfg.domain == "ftir"
                  else r"2$\theta$ / deg")
    ax.set_ylabel(f"{tag} loading (a.u.)")
    ax.legend(fontsize=6.5)
    return ax


# ===================================================== the composite deliverable
def diagnostics_figure(X, y, x, cfg, groups=None, pre_options=None, ref=None,
                       ref_label="reference", highlight=None, viz_pre=None,
                       viz_components=None, basename="pls_diagnostics", method="pls"):
    """The 2x2 PLS diagnostic the family exists to produce, gated and QA'd:
       (a) RMSECV vs components, one curve per preprocessing, parsimony pick ringed
       (b) regression-coefficient spectrum (+ optional pure-component reference)
       (c) scores coloured by concentration
       (d) loadings of the first two components

    The honest error (a) is computed from a leakage-safe scan; (b)-(d) are drawn from
    a model fit on the whole set for interpretability. By default that model uses the
    parsimony pick, but pass viz_pre / viz_components to show a more band-interpretable
    model (e.g. raw, 2 LV) while still REPORTING the SNV/1-LV error in (a).

    Returns (fig, axes, info) where info has the scan, the parsimony choice, and the
    fitted viz model. Pass groups=level-labels to make (a) a leave-one-level-out error."""
    import matplotlib.pyplot as plt
    _, X, y = _check_xy(X, y, cfg)
    x = np.asarray(x, float)
    if pre_options is None:
        pre_options = tuple(dict.fromkeys((cfg.pls_preprocess, "none")))
    scan = component_scan(X, y, cfg, pre_options=pre_options, groups=groups, method=method)
    choice = choose_n_components(scan, cfg)
    fitter = {"pls": pls_fit, "pcr": pcr_fit}[method]
    model = fitter(X, y, viz_components or choice["n_components"], cfg,
                   pre=viz_pre or choice["pre"])

    style.apply_style(cfg)
    fig, ax = style.figure(cfg, 2, 2, layout="constrained")   # colorbar must be layout-aware
    (a, b), (c, d) = ax
    plot_rmsecv(a, scan, cfg, choice=choice)
    a.set_title("How many components?", fontsize=8.5)
    plot_coefficients(b, model, x, cfg, ref=ref, ref_label=ref_label, highlight=highlight)
    b.set_title(f"Coefficients ({model['pre_name']})", fontsize=8.5)
    plot_scores(c, model, cfg, color_by=y, cbar=fig)
    c.set_title("Scores", fontsize=8.5)
    comps = (0, 1) if model["scores"].shape[1] > 1 else (0,)
    plot_loadings(d, model, x, cfg, comps=comps, highlight=highlight)
    d.set_title("Loadings", fontsize=8.5)

    style.finalize_figure(fig, wspace=0.10, hspace=0.12)       # constrained-layout fractions
    # left column clears its own narrow y-furniture; right column gets a tight offset so
    # b/d sit at their corners instead of being dragged left by d's wide decimal ticks.
    style.add_panel_labels(fig, cfg, axes=[a, c], labels=["a", "c"], x_offset_pt="auto")
    style.add_panel_labels(fig, cfg, axes=[b, d], labels=["b", "d"], x_offset_pt=-10.0)
    import os
    verify.render_preview(fig, os.path.join(cfg.output_dir, "_preview_" + basename + ".png"))
    verify.audit_layout(fig, cfg)
    return fig, (a, b, c, d), {"scan": scan, "choice": choice, "model": model}


# ===================================================== methods text
def methods_text(model, cv, choice, cfg):
    """Paste-ready methods paragraph: matrix, preprocessing, model + component count
    (with the selection rule), CV scheme + the question it answers, and the figures of
    merit. Edit to taste; do not invent numbers it doesn't contain."""
    n = model["y"].size
    p = np.asarray(model["coef"]).size
    tag = "latent variables" if model["kind"] == "pls" else "principal components"
    lines = [
        f"Multivariate calibration was performed by {model['kind'].upper()} on the "
        f"baseline-corrected spectra ({n} samples x {p} wavenumbers).",
        f"Spectra were preprocessed by {model['pre_name']} before modelling.",
        f"The number of {tag} was chosen by cross-validation: {choice['caption']}.",
        f"Cross-validation used {cv['scheme']}, which answers '{cv['question']}'.",
        f"Figures of merit: RMSEC={model.get('rmsec', float('nan')):.3g} "
        f"{cfg.conc_unit}, R2(cal)={model.get('r2_cal', float('nan')):.3f}; "
        f"RMSECV={cv['rmsecv']:.3g} {cfg.conc_unit}, "
        f"R2(CV)={cv['r2cv']:.3f}, bias={cv['bias']:+.3g}.",
    ]
    return "\n".join(lines)


# ===================================================== outlier / diagnostic gate
# PCA-based per-sample diagnostics, the modern chemometric outlier suite:
#   * Hotelling T2  (score distance) -- "extreme but maybe valid" position in the model
#   * Q / SPE       (orthogonal distance) -- "doesn't fit the model" (the usual bad-load tell)
#   * leverage      (Sigma U^2; the unknown-sample extrapolation flag at 2A/n or 3A/n)
#   * DD-SIMCA      (Pomerantsev & Rodionova 2014) -- SD + OD combined into a chi^2(Nh+Nq)
#                   statistic with an extreme limit (alpha) and a Bonferroni outlier limit (gamma)
# numpy + scipy only (NO sklearn) so it runs as raw-spectrum QC BEFORE any calibration. The
# limit forms are exactly those in the literature (Jackson & Mudholkar 1979 for Q; the F / chi^2
# / Tracy-Young-Mason 1992 forms for T2). Validated bit-for-bit against the hand-derived fixture
# in skill_validation/chemometrics/golden/closed_form.py, plus scipy F/chi2 closed forms.
#
# Why a separate PCA decomposition rather than reusing pca_fit's dict: this is intentionally
# self-contained (numpy SVD), stores the training mean/std so NEW unknowns can be scored, and
# matches the closed-form derivation exactly. T2 = (n-1)*leverage; mean(leverage) = A/n.

def pca_outlier_model(X, cfg=None, n_components=2, pre="none", scale=False):
    """Decompose X for outlier diagnostics: preprocess (`pre`), centre (store the mean),
    optionally autoscale (ddof=1), SVD, retain A = min(n_components, numerical rank).
    `pre` defaults to "none" (raw-spectrum QC); pass e.g. "snv" to QC the modelled space.
    Stores the training mean/std so `outlier_distances(om, Xnew)` can score new samples."""
    X = np.asarray(X, float)
    pp = Preprocessor(pre, cfg)
    Xp = pp.fit_transform(X)
    mean = Xp.mean(0)
    Xc = Xp - mean
    std = None
    if scale:
        sd = Xc.std(0, ddof=1)
        std = np.where(sd == 0, 1.0, sd)
        Xc = Xc / std
    n, p = Xc.shape
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    tol = (s.max() * max(n, p) * np.finfo(float).eps) if s.size else 0.0
    rank = int(np.sum(s > tol))
    A = int(max(1, min(n_components, rank)))
    eig = (s ** 2) / (n - 1)                        # eigenvalues (ddof=1)
    return {"kind": "pca_outlier", "pre": pp, "mean": mean, "std": std, "scale": scale,
            "s": s, "eig": eig, "P": Vt.T, "U": U, "scores": U * s,
            "A": A, "rank": rank, "n": n, "p": p}


def _outlier_centred(om, X):
    if X is None:
        return om["scores"] @ om["P"].T            # exact reconstruction of the training Xc
    Xp = om["pre"].transform(np.asarray(X, float))
    Xc = Xp - om["mean"]
    if om["scale"] and om["std"] is not None:
        Xc = Xc / om["std"]
    return Xc


def outlier_distances(om, X=None):
    """Per-sample T2, Q/SPE, leverage and the DD-SIMCA score/orthogonal distances, for the
    training set (X=None) or for new samples X. Q = ||Xc||^2 - sum(T_A^2) is general (it also
    captures the part of a new sample outside the model span). leverage = sum (T_a/s_a)^2."""
    A = om["A"]; n = om["n"]
    s = om["s"][:A]
    eig = om["eig"][:A]
    Xc = _outlier_centred(om, X)
    T = Xc @ om["P"]
    Ta = T[:, :A]
    t2 = np.sum(Ta ** 2 / np.where(eig == 0, np.inf, eig), axis=1)
    q = np.sum(Xc ** 2, axis=1) - np.sum(Ta ** 2, axis=1)
    q = np.clip(q, 0.0, None)                       # numerical floor (perfect fit -> 0)
    lev = np.sum((Ta / np.where(s == 0, np.inf, s)) ** 2, axis=1)   # = T2/(n-1); mean = A/n
    return {"t2": t2, "q": q, "leverage": lev, "sd": lev, "od": q, "scores": T}


def t2_limit(om, alpha=0.05, form="F_textbook"):
    """Hotelling T2 critical limit. forms: 'F_textbook' [A(n-1)/(n-A)]F (mdatools),
    'F_n2m1' [A(n^2-1)/(n(n-A))]F (chemometrics lib), 'chi2' (large-n), 'beta'
    (Tracy-Young-Mason, Phase-I). Falls back to chi2 when n-A<=0."""
    from scipy.stats import f as fdist, chi2, beta
    A = om["A"]; n = om["n"]
    if form == "chi2" or n - A <= 0:
        return float(chi2.ppf(1 - alpha, A))
    fq = float(fdist.ppf(1 - alpha, A, n - A))
    if form == "F_textbook":
        return A * (n - 1) / (n - A) * fq
    if form == "F_n2m1":
        return A * (n ** 2 - 1) / (n * (n - A)) * fq
    if form == "beta":
        if n - A - 1 <= 0:
            return float(chi2.ppf(1 - alpha, A))
        return (n - 1) ** 2 / n * float(beta.ppf(1 - alpha, A / 2.0, (n - A - 1) / 2.0))
    raise ValueError(f"unknown T2 limit form: {form!r}")


def q_limit(om, alpha=0.05, method="jackson_mudholkar", q_values=None):
    """Q/SPE critical limit. 'jackson_mudholkar' uses the discarded eigenvalues (theta1-3,h0);
    'box' fits g*chi2_h to the moments of the Q values (needs q_values). Returns 0.0 when there
    is no residual space (A>=rank) or a perfect fit -- the Q test is then uninformative."""
    from scipy.stats import norm, chi2
    if method == "jackson_mudholkar":
        res = om["eig"][om["A"]:om["rank"]]
        res = res[res > 0]
        if res.size == 0:
            return 0.0
        th1 = float(res.sum()); th2 = float((res ** 2).sum()); th3 = float((res ** 3).sum())
        if th2 == 0:
            return 0.0
        h0 = 1.0 - 2.0 * th1 * th3 / (3.0 * th2 ** 2)
        if h0 <= 0:
            h0 = 1e-6
        z = float(norm.ppf(1 - alpha))
        return float(th1 * (z * np.sqrt(2.0 * th2 * h0 ** 2) / th1 + 1.0
                            + th2 * h0 * (h0 - 1.0) / th1 ** 2) ** (1.0 / h0))
    if method == "box":
        if q_values is None:
            raise ValueError("box method needs q_values")
        qv = np.asarray(q_values, float)
        m = float(qv.mean()); v = float(qv.var(ddof=1))
        if m <= 0 or v <= 0:
            return 0.0
        g = v / (2.0 * m); h = 2.0 * m ** 2 / v
        return float(g * chi2.ppf(1 - alpha, h))
    raise ValueError(f"unknown Q limit method: {method!r}")


def _dd_dof(d, method="classical"):
    """DD-SIMCA data-driven degrees of freedom + scaling for a distance vector, under the
    parameterization u ~ (u0/N)*chi2(N) (so N*u/u0 ~ chi2(N)).

    method='classical' (mdatools ddmoments): N = round(2*mean^2/var), u0 = mean. EXACT but
        NON-ROBUST — a strong outlier inflates mean/var and masks itself. Right for clean data.
    method='robust' (mdatools ddrobust; Pomerantsev & Rodionova 2014): solve N from the
        IQR/median ratio of chi2(N) (which is scale-free), then u0 = N*median/median(chi2_N).
        Use this for DETECTION — the median/IQR are not dragged by the outliers being sought.
    Both clamped to [1,250]."""
    from scipy.stats import chi2
    d = np.asarray(d, float)
    if method == "robust":
        M = float(np.median(d))
        R = float(np.percentile(d, 75) - np.percentile(d, 25))
        if M <= 0 or R <= 0:
            return 1, max(M, 1e-12)
        ratio = R / M

        def g(nu):
            return (chi2.ppf(0.75, nu) - chi2.ppf(0.25, nu)) / chi2.ppf(0.5, nu) - ratio
        lo, hi = 1.0, 250.0
        glo, ghi = g(lo), g(hi)
        if glo * ghi < 0:
            try:
                from scipy.optimize import brentq
                nu = brentq(g, lo, hi, xtol=1e-4)
            except Exception:
                nu = lo
        else:
            nu = lo if abs(glo) <= abs(ghi) else hi      # ratio outside the chi2 range -> clamp
        N = max(1, min(250, int(round(nu))))
        u0 = N * M / float(chi2.ppf(0.5, N))
        return N, u0
    m = float(d.mean()); v = float(d.var(ddof=1))
    N = int(round(2.0 * m ** 2 / v)) if v > 0 else 1
    return max(1, min(250, N)), m


def ddsimca_limits(om, sd, od, alpha=0.05, gamma=0.01, dof="classical"):
    """DD-SIMCA combined-distance limits. f = Nh*(SD/h0) + Nq*(OD/q0) ~ chi2(Nh+Nq); extreme
    limit at alpha; per-sample Bonferroni outlier limit at gamma (None => no outlier line).
    `dof`: 'classical' (moments) or 'robust' (median/IQR — use for detection, see _dd_dof)."""
    from scipy.stats import chi2
    n = om["n"]
    Nh, h0 = _dd_dof(sd, dof)
    od = np.asarray(od, float)
    if od.size == 0 or float(np.max(od)) <= 1e-12 * max(1.0, float(np.mean(sd))):
        print("  [WARN] ddsimca: residual space is empty (components >= rank of the training set) - "
              "the orthogonal distance carries no information; limits use the score distance only")
        Nq, q0 = 0, float("inf")                       # Nq*(OD/q0) vanishes instead of 0/0
    else:
        Nq, q0 = _dd_dof(od, dof)
    c_crit = float(chi2.ppf(1 - alpha, Nh + Nq))
    c_out = float(chi2.ppf((1 - gamma) ** (1.0 / n), Nh + Nq)) if gamma else None
    return {"Nh": Nh, "Nq": Nq, "h0": h0, "q0": q0, "c_crit": c_crit, "c_out": c_out, "dof": dof}


def diagnostics(X, cfg=None, n_components=2, pre="none", scale=False, alpha=0.05, gamma=0.01,
                t2_form="F_textbook", q_method="jackson_mudholkar", lev_factor=2, dof="robust"):
    """The full outlier/diagnostic read-out for a spectral matrix: builds the PCA model,
    computes T2/Q/leverage + their limits, the DD-SIMCA combined distance + flags
    (regular/extreme/outlier). Returns one dict; feed it to plot_influence. The DD-SIMCA
    path is the primary verdict; T2(F) and Q(Jackson-Mudholkar) are the classical cross-checks.
    `dof` defaults to 'robust' (median/IQR) so strong outliers don't mask themselves; pass
    'classical' for the method-of-moments DoF (right only when the set is already clean)."""
    om = pca_outlier_model(X, cfg, n_components, pre, scale)
    d = outlier_distances(om)
    dd = ddsimca_limits(om, d["sd"], d["od"], alpha, gamma, dof=dof)
    f = dd["Nh"] * (d["sd"] / dd["h0"]) + dd["Nq"] * (d["od"] / dd["q0"])
    if dd["c_out"] is not None:
        flags = np.where(f > dd["c_out"], "outlier",
                         np.where(f > dd["c_crit"], "extreme", "regular"))
    else:
        flags = np.where(f > dd["c_crit"], "extreme", "regular")
    return {"model": om, **d,
            "t2_limit": t2_limit(om, alpha, t2_form), "t2_form": t2_form,
            "q_limit": q_limit(om, alpha, q_method, q_values=d["q"]), "q_method": q_method,
            "leverage_warn": lev_factor * om["A"] / om["n"],
            "dd": dd, "f": f, "flags": flags, "alpha": alpha, "gamma": gamma,
            "caption": (f"PCA outlier diagnostics ({om['A']} PC, {om['pre'].name}): "
                        f"DD-SIMCA[{dd['dof']}] chi^2({dd['Nh']}+{dd['Nq']}), extreme alpha={alpha}"
                        + (f", outlier gamma={gamma}" if gamma else "")
                        + f"; {int(np.sum(flags!='regular'))} of {om['n']} flagged.")}


def plot_influence(ax, diag, cfg=None, labels=None, axis_scale="linear"):
    """DD-SIMCA acceptance ('influence') plot: orthogonal distance OD/q0 vs score distance
    SD/h0, with the chi^2 boundary drawn as a straight line (extreme dashed, outlier dotted).
    Points keyed regular/extreme/outlier by marker+colour (grayscale-safe). Annotate by
    passing labels (e.g. sample ids) — only flagged points get a tag.

    axis_scale: 'linear' (default) shows the true linear boundary — best when the distances are
    of comparable size. Use 'sqrt' when one or a few GROSS outliers dominate the range and squash
    the in-control cluster and the limit lines into the corner (the single-bad-spectrum case): a
    sqrt axis compresses the large values so the borderline/extreme points and the boundary stay
    legible. The boundary stays correct — it is drawn in data coordinates and the axis transform
    curves it; sqrt (not log) is the default transform because SD/OD are >=0 and often exactly 0."""
    dd = diag["dd"]
    xs = np.asarray(diag["leverage"]) / dd["h0"]
    ys = np.asarray(diag["q"]) / dd["q0"]
    flags = np.asarray(diag["flags"])
    sty = {"regular": ("#0072B2", "o"), "extreme": ("#E69F00", "s"), "outlier": ("#D55E00", "D")}
    for fl, (col, mk) in sty.items():
        m = flags == fl
        if np.any(m):
            ax.scatter(xs[m], ys[m], c=col, marker=mk, s=30, edgecolor="white",
                       linewidth=0.3, label=fl, zorder=3)
    xmax = max(float(xs.max()) * 1.1, dd["c_crit"] / dd["Nh"]) if xs.size else 1.0
    xl = np.linspace(0, xmax, 100)
    for c, ls, lab in [(dd["c_crit"], "--", f"extreme (α={diag['alpha']})"),
                       (dd["c_out"], ":", f"outlier (γ={diag['gamma']})" if dd["c_out"] else None)]:
        if c is None:
            continue
        yb = (c - dd["Nh"] * xl) / dd["Nq"]
        ok = yb >= 0
        ax.plot(xl[ok], yb[ok], ls, color="0.35", lw=0.9, label=lab, zorder=2)
    if labels is not None:
        for i, fl in enumerate(flags):
            if fl != "regular":
                ax.annotate(str(labels[i]), (xs[i], ys[i]), fontsize=6,
                            xytext=(2, 2), textcoords="offset points")
    if axis_scale == "sqrt":
        fwd = lambda a: np.sqrt(np.clip(np.asarray(a, float), 0, None))
        inv = lambda a: np.asarray(a, float) ** 2
        ax.set_xscale("function", functions=(fwd, inv))
        ax.set_yscale("function", functions=(fwd, inv))
    elif axis_scale not in ("linear", None):
        raise ValueError(f"axis_scale must be 'linear' or 'sqrt', got {axis_scale!r}")
    ax.set_xlabel("score distance / $h_0$ (ratio)")
    ax.set_ylabel("orthogonal distance / $q_0$ (ratio)")
    ax.legend(fontsize=6, loc="best")
    return ax


# =============================================== DD-SIMCA one-class classifier + class FoM
# The outlier gate above asks "is this sample weird for THIS set?"; a one-class CLASSIFIER asks
# "is this sample a member of the TARGET class?" — the right frame for cocrystal ID ("is this the
# target phase, or a physical mixture / a different polymorph?"). Fit the DD-SIMCA acceptance model
# on the target class's spectra, then accept/reject unknowns against its chi^2 boundary, and report
# the class figures-of-merit (sensitivity / specificity / efficiency) — NOT R2/LOD, the wrong
# vocabulary for a classifier. Reuses the validated pca_outlier_model / ddsimca_limits path, so its
# limits inherit test_outliers' validation. Routed to from references/cocrystal_id.md.

def ddsimca_fit(X, cfg=None, n_components=2, pre="none", scale=False,
                alpha=0.05, gamma=None, dof="classical"):
    """Fit a DD-SIMCA ONE-CLASS model on the TARGET class's spectra: a PCA model + the SD/OD
    scaling (h0, q0, Nh, Nq) + the acceptance limit c_crit at 1-alpha. Feed the model to
    `ddsimca_predict` to accept/reject unknowns ('is this the target phase?'). `dof` defaults to
    'classical' — a training class is assumed clean; use 'robust' only if it may itself contain
    outliers. gamma=None → no per-sample outlier line (a class model needs only the acceptance
    boundary)."""
    om = pca_outlier_model(X, cfg, n_components, pre, scale)
    d = outlier_distances(om)
    dd = ddsimca_limits(om, d["sd"], d["od"], alpha, gamma, dof=dof)
    return {"kind": "ddsimca", "om": om, "dd": dd, "alpha": alpha,
            "train_sd": d["sd"], "train_od": d["od"],
            "caption": (f"DD-SIMCA one-class model ({om['A']} PC, {om['pre'].name}): accept if the "
                        f"chi^2({dd['Nh']}+{dd['Nq']}) combined distance <= c_crit (alpha={alpha}).")}


def ddsimca_predict(model, X):
    """Accept/reject samples X against a `ddsimca_fit` model. Returns dict: `f` (combined
    distance), `accept` (bool — f <= c_crit → IN the target class), `sd`, `od`, `c_crit`. A
    rejected unknown is NOT the modelled class: a target-cocrystal model rejects a physical
    mixture or a different polymorph. X is projected through the STORED model (out-of-sample)."""
    om = model["om"]; dd = model["dd"]
    dn = outlier_distances(om, np.atleast_2d(np.asarray(X, float)))
    f = dd["Nh"] * (dn["sd"] / dd["h0"]) + dd["Nq"] * (dn["od"] / dd["q0"])
    return {"f": f, "accept": f <= dd["c_crit"], "sd": dn["sd"], "od": dn["od"],
            "c_crit": dd["c_crit"]}


def class_fom(y_true, y_pred):
    """Classification figures-of-merit — the CORRECT vocabulary for a class model (DD-SIMCA /
    SIMCA / PLS-DA), not R2/LOD. `y_true`, `y_pred` are boolean arrays (True = member of the
    target class / accepted):
      sensitivity (TPR) = TP/(TP+FN)  — target members correctly accepted,
      specificity (TNR) = TN/(TN+FP)  — non-members correctly rejected,
      efficiency        = sqrt(sensitivity·specificity)  — the SIMCA class efficiency.
    Returns those + the confusion counts + a caption. NaN-safe (a rate with no denominator is nan)."""
    yt = np.asarray(y_true, bool).ravel()
    yp = np.asarray(y_pred, bool).ravel()
    if yt.shape != yp.shape:
        raise ValueError("y_true and y_pred must have the same shape")
    tp = int(np.sum(yt & yp)); fn = int(np.sum(yt & ~yp))
    tn = int(np.sum(~yt & ~yp)); fp = int(np.sum(~yt & yp))
    sens = tp / (tp + fn) if (tp + fn) else float("nan")
    spec = tn / (tn + fp) if (tn + fp) else float("nan")
    eff = float(np.sqrt(sens * spec)) if (np.isfinite(sens) and np.isfinite(spec)) else float("nan")
    caption = (f"class FoM — sensitivity={sens:.3f} (TPR {tp}/{tp+fn}), "
               f"specificity={spec:.3f} (TNR {tn}/{tn+fp}), efficiency={eff:.3f} "
               f"(geom. mean) [SIMCA class metrics, not R2/LOD]")
    return {"sensitivity": sens, "specificity": spec, "efficiency": eff,
            "tp": tp, "fn": fn, "tn": tn, "fp": fp, "caption": caption}


# ===================================================== validation trio
# Cheap, defensible figures-of-merit for a SMALL calibration: a permutation/y-scramble test
# (is the model better than chance?), a corrected paired t-test (is model A really better than
# B?), and RPD/RPIQ (is the calibration useful?).

def _permute_y(y, groups, rng):
    """Shuffle y for the null. With groups (level labels) permute BETWEEN levels (each level's
    reps keep ONE shared y) — replicate structure is preserved; sklearn's groups= permutes
    WITHIN groups, which is wrong for replicate spectra. Without groups, a plain full shuffle."""
    y = np.asarray(y, float)
    if groups is None:
        return rng.permutation(y)
    groups = np.asarray(groups)
    uniq = list(dict.fromkeys(groups.tolist()))
    vals = rng.permutation(np.array([y[groups == g].mean() for g in uniq]))
    yp = np.empty_like(y)
    for g, v in zip(uniq, vals):
        yp[groups == g] = v
    return yp


def permutation_test(X, y, n_components, cfg, n_perm=199, groups=None, scheme=None,
                     method="pls", pre=None, seed=0):
    """Permutation / y-scrambling test: refit on shuffled y n_perm times to build the null for
    R2(CV), then p = (#{null >= observed} + 1) / (n_perm + 1). With groups, y is permuted between
    LEVELS (see _permute_y). Guards a small calibration against a chance correlation.

    The per-fold preprocessing depends on X only, so it is done once and reused across the
    permutations whenever the folds themselves do not depend on y (groups given, or the
    per-sample LOO scheme). The level-stratified k-fold WITHOUT groups builds its folds from
    the y values, so there the folds and blocks are rebuilt per permutation."""
    rng = np.random.default_rng(seed)
    pre = cfg.pls_preprocess if pre is None else pre
    _, X, y = _check_xy(X, y, cfg, n_components)
    obs = cross_validate(X, y, n_components, cfg, pre=pre, groups=groups, scheme=scheme,
                         method=method, report=False)["r2cv"]
    folds_fixed = groups is not None or _scheme(scheme, cfg) in ("auto", "loo")
    if folds_fixed:
        fixed_folds, _, _ = _make_folds(len(y), y, groups, scheme, cfg)
        fixed_blocks = _fold_blocks(X, cfg, pre, fixed_folds)
    null = np.empty(n_perm)
    for i in range(n_perm):
        yp = _permute_y(y, groups, rng)
        folds = fixed_folds if folds_fixed else _make_folds(len(yp), yp, groups, scheme, cfg)[0]
        blocks = fixed_blocks if folds_fixed else None
        pred = _cv_predict(X, yp, n_components, cfg, pre, folds, method, blocks=blocks)
        null[i] = _cv_stats(yp, pred)[2]
    p = (int(np.sum(null >= obs)) + 1) / (n_perm + 1)
    return {"observed": float(obs), "null": null, "p_value": float(p), "n_perm": n_perm,
            "caption": f"permutation test: R2(CV)={obs:.3f}, p={p:.3g} ({n_perm} permutations)"}


def corrected_paired_t(scores_a, scores_b, n_train, n_test):
    """Nadeau & Bengio (2003) CORRECTED resampled paired t-test comparing two models' per-resample
    scores/errors. The naive paired-t variance is inflated by (1/k + n_test/n_train) to account
    for the training sets overlapping across resamples. Returns t, two-sided p, df=k-1, correction.
    (mlxtend's paired_ttest_resampled is the UNcorrected Dietterich formula — not this.)"""
    from scipy.stats import t as tdist
    a = np.asarray(scores_a, float); b = np.asarray(scores_b, float)
    d = a - b
    k = d.size
    md = float(d.mean()); var = float(d.var(ddof=1))
    correction = 1.0 / k + float(n_test) / float(n_train)
    se = np.sqrt(correction * var)
    if se > 0:
        t = md / se
    else:
        t = 0.0 if md == 0 else float(np.inf) * np.sign(md)
    p = float(2 * tdist.sf(abs(t), k - 1)) if np.isfinite(t) else 0.0
    return {"t": float(t), "p_value": p, "df": k - 1, "correction": float(correction),
            "mean_diff": md,
            "caption": f"corrected paired t = {t:.3f}, p = {p:.3g} (df={k-1}, "
                       f"correction={correction:.3f}; Nadeau-Bengio)"}


def rpd(y_reference, rmse, ddof=1):
    """Ratio of Performance to Deviation = SD(reference)/RMSE. SD uses ddof=1 (sample). Rough
    interpretation (Chang 2001): >2 useful, >2.5 good, >3 excellent. Pass RMSEP or RMSECV."""
    sd = float(np.std(np.asarray(y_reference, float), ddof=ddof))
    return float(sd / rmse) if rmse > 0 else float("inf")


def rpiq(y_reference, rmse):
    """Ratio of Performance to IQR = IQR(reference)/RMSE (IQR via R-type-7 / numpy default
    linear interpolation). RPIQ ~= 1.349*RPD for normal y; more robust to skew."""
    y = np.asarray(y_reference, float)
    iqr = float(np.percentile(y, 75) - np.percentile(y, 25))
    return float(iqr / rmse) if rmse > 0 else float("inf")
