"""The dependency-free closed-form outlier-diagnostics fixture (SPEC.md §7a).

A 6x3 matrix, centre-only, A=2 components, alpha=0.05. Every expected value was
hand-derived and independently reproduced in scipy during the research pass. This is the
always-runs anchor for the outlier gate (step 2): it needs no dataset and no sklearn.
The DD-SIMCA moment DoF at n=6 are noisy and are RECOMPUTED + locked when the gate is
implemented (left as None here).
"""
import numpy as np

X = np.array([
    [2.0, 0.0, 1.0],
    [0.0, 2.0, 0.0],
    [3.0, 1.0, 2.0],
    [1.0, 3.0, 1.0],
    [4.0, 0.0, 3.0],
    [2.0, 2.0, 2.0],
])
N, M = X.shape          # 6, 3
A = 2                   # retained components
ALPHA = 0.05
PREPROCESS = "center only (no scaling)"

EIGENVALUES_DDOF1 = [3.7812886, 0.7632012, 0.0221769]

T2_PER_SAMPLE = [2.4705981, 2.2419840, 0.3937636, 1.9042242, 2.1356383, 0.8537918]
T2_LIMITS = {
    "F_n2m1": 20.2541264,    # [k(n^2-1)/(n(n-k))] * F   (chemometrics lib form)
    "F_textbook": 17.3606798,  # [k(n-1)/(n-k)] * F        (mdatools form)
    "chi2": 5.9914645,
    "beta_tracy_young": 3.6011630,
}

Q_PER_SAMPLE = [0.00065204, 0.00572200, 0.04670982, 0.01321247, 0.00808042, 0.03650782]
Q_LIMITS = {
    "jackson_mudholkar": 0.08309166,   # degenerate single residual eigenvalue, h0=1/3
    "box_nomikos_macgregor": 0.05571389,   # g=0.009407, h=1.964643
}

DDSIMCA = {
    "mean_SD_equals_A_over_N": 0.3333,   # confirms singular-value normalization of SD
    # locked from the implementation (method-of-moments DoF, ddof=1), 2026-06-23:
    "Nh": 8, "Nq": 2,
    "h0": 0.3333333333333333,             # mean score distance (= A/n)
    "q0": 0.018480761891874014,           # mean orthogonal distance
    "Dcrit": 18.307038053275146,          # chi2(0.95, Nh+Nq=10)  -- extreme limit (alpha=0.05)
    "c_out": 28.20507535399202,           # chi2((1-0.01)^(1/6), 10) -- Bonferroni outlier (gamma=0.01)
}
