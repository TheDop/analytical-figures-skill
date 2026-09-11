"""pytest config for the skill self-tests: headless matplotlib + skill root on sys.path."""
import os
import sys
import matplotlib

matplotlib.use("Agg")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # the skill root
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
