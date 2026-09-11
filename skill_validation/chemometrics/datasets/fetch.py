"""Optional dataset fetcher for the cross-language golden tests.

Tests that need an external dataset SKIP cleanly when it is absent, so this is only
needed to upgrade those SKIPs to real T2 checks. Run it from a machine with network access.

  gasoline  -> datasets/gasoline.npz  (X: 60x401 NIR, y: octane) -- R `pls` benchmark
  emsc      -> datasets/emsc_testdata.xlsx                       -- biospectools EMSC gold

We do not vendor these by default (licensing / size). Sources:
  gasoline: R `pls` package `data(gasoline)` (Kalivas 1997). Export from R:
      library(pls); data(gasoline)
      write.csv(cbind(octane=gasoline$octane, gasoline$NIR), "gasoline.csv", row.names=FALSE)
    then load the CSV here and save as gasoline.npz.
  emsc: https://raw.githubusercontent.com/BioSpecNorway/biospectools/main/tests/data/emsc_testdata.xlsx
"""
from __future__ import annotations
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def fetch_emsc():
    import urllib.request
    url = ("https://raw.githubusercontent.com/BioSpecNorway/biospectools/"
           "main/tests/data/emsc_testdata.xlsx")
    dst = HERE / "emsc_testdata.xlsx"
    urllib.request.urlretrieve(url, dst)
    print(f"wrote {dst}")


def gasoline_from_csv(csv_path):
    """Convert an exported gasoline.csv (first column octane, rest NIR) -> gasoline.npz."""
    import numpy as np
    raw = np.genfromtxt(csv_path, delimiter=",", skip_header=1)
    y, X = raw[:, 0], raw[:, 1:]
    np.savez(HERE / "gasoline.npz", X=X, y=y)
    print(f"wrote {HERE / 'gasoline.npz'}  X{X.shape} y{y.shape}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "emsc":
        fetch_emsc()
    elif len(sys.argv) > 2 and sys.argv[1] == "gasoline":
        gasoline_from_csv(sys.argv[2])
    else:
        print(__doc__)
