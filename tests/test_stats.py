"""Multiplicity correction matches statsmodels (Bonferroni / Holm / BH-FDR)."""
import numpy as np
import pytest
from scripts import verify

PVALS = np.array([0.001, 0.013, 0.04, 0.2, 0.5, 0.9])


@pytest.mark.parametrize("method,sm_method", [("bonferroni", "bonferroni"),
                                              ("holm", "holm"), ("bh", "fdr_bh")])
def test_adjust_vs_statsmodels(method, sm_method):
    mt = pytest.importorskip("statsmodels.stats.multitest")
    got = verify.adjust_pvalues(PVALS, method)
    want = mt.multipletests(PVALS, method=sm_method)[1]
    assert np.allclose(got, want, atol=1e-12), f"{method}: {got} vs {want}"


def test_order_preserved():
    shuffled = PVALS[::-1]
    adj = verify.adjust_pvalues(shuffled, "holm")
    assert adj.shape == shuffled.shape
    # adjusting then re-sorting equals adjusting the sorted vector
    assert np.allclose(np.sort(adj), np.sort(verify.adjust_pvalues(PVALS, "holm")))


def test_multiplicity_check_flags_raw_inflation():
    # several raw-significant p-values; correction should drop the count -> result is coherent
    res = verify.multiplicity_check([0.001, 0.02, 0.03, 0.045], method="holm", alpha=0.05)
    assert res["adjusted"].max() <= 1.0
    assert int(res["reject"].sum()) <= int(np.sum(np.asarray(res["raw"]) <= 0.05))
