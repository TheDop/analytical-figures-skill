"""Golden reference numbers captured as LITERALS (so the suite is R-free at test time).

Each entry: value(s) + source + version/settings + tier. Used by test_* modules. Numbers
sourced from package vignettes / committed unit tests / papers during the 2026-06-23
research pass; see SPEC.md §7 and §10 for provenance. Entries marked NEW are consumed once
the corresponding feature (outlier gate / EMSC / trio) is implemented (step 2).
"""

# --- R `pls` gasoline, summary(plsr(octane~NIR, ncomp=10, validation="LOO")) CV row [T2]
GASOLINE_RMSEP_CV = {  # ncomp -> RMSEP (LOO)
    0: 1.545, 1: 1.357, 2: 0.2966, 3: 0.2524, 4: 0.2476, 5: 0.2398,
    6: 0.2319, 7: 0.2386, 8: 0.2316, 9: 0.2449, 10: 0.2673,
}
GASOLINE_PREDICT_2C = [87.94, 87.25, 88.16, 84.97, 85.15, 84.51, 87.56, 86.85, 89.19, 87.09]
GASOLINE_SOURCE = "R pls vignette; data=gasoline (Kalivas 1997); train[1:50]/test[51:60]; scale=FALSE"

# --- mdatools people, pca(people, 11, scale=TRUE) Hotelling T2 extreme limits per comp [T2, NEW]
MDATOOLS_PEOPLE_T2LIM = [4.15961510, 6.85271430, 9.40913034]   # comps 1..3 (comp 11 = 37.07020869)
MDATOOLS_PEOPLE_DDMOMENTS_C1 = dict(Nq=10, Nh=2, q0=5.39623603, h0=0.96875)
MDATOOLS_SOURCE = "mdatools tests/testthat/test-pca.R; data(people); scale=TRUE; alpha=0.05/gamma=0.01"

# --- chemometrics (Rüdt) peach NIR, PLSRegression(n_components=4), center only [T1, NEW]
CHEMOMETRICS_PEACH = dict(x_residual_std=0.005448557285711451,
                          crit_dmodx_095=1.0697876760,
                          crit_dhypx_095=11.1869800486,
                          sum_leverage=4.0)
CHEMOMETRICS_SOURCE = "maruedt/chemometrics v0.4.0; bundled peach Brix NIR (50x600)"

# --- PyChemAuth DD-SIMCA, iris-setosa (25x4), n_components=3, scale_x=False, robust='semi' [T2, NEW]
PYCHEMAUTH_SETOSA_DDSIMCA = dict(Nh=3, Nq=1, h0=3.2926643176398445, q0=0.005740692404315845,
                                 c_crit=9.487729036781154)   # chi2_0.95,4
PYCHEMAUTH_SOURCE = "mahynski/pychemauth tests/test_ddsimca.py; benchmarks mdatools 0.14.1"

# --- EMSC biospectools gold (poly_order=4, reference=mean) [T1, NEW]
EMSC_BIOSPECTOOLS = dict(file="datasets/emsc_testdata.xlsx", poly_order=4, reference="mean",
                         max_abs_err_target=1e-12)   # parity to MATLAB was 6e-16
EMSC_SOURCE = ("BioSpecNorway/biospectools tests/data/emsc_testdata.xlsx "
               "(validated vs Kohler MATLAB 06.02.2020); Afseth & Kohler 2012")

# --- validation trio goldens [T0, NEW]
NADEAU_BENGIO_GOLDEN = dict(t=1.364576, p=0.205528, df=9, correction=0.35)
NADEAU_BENGIO_SOURCE = "RNG-free golden vector; Nadeau & Bengio 2003; cross-checked vs R correctR"
PERMUTATION_GOLDEN = dict(score=0.829, p=0.00498, n_perm=200)   # recipe in SPEC; p=(#>=obs+1)/(n_perm+1)
PERMUTATION_SOURCE = "sklearn permutation_test_score mechanics; p=(#{perm>=obs}+1)/(n_perm+1)"
RPD_UNIT = dict(sd=10.0, rmsep=4.0, rpd=2.50)   # RPD = SD(ddof1)/RMSEP
RPD_SOURCE = "Williams; Bellon-Maurel 2010 (RPIQ~1.349*RPD); Chang 2001 thresholds"
