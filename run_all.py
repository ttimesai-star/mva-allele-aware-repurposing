"""Reproduce every output of this repository: python run_all.py [--offline]

--offline uses only the committed HTTP cache in data/raw/http_cache (no network)."""
import subprocess
import sys

STEPS = ["verify_refs.py", "clinvar_landscape.py", "drug_annotations.py", "score.py", "make_figures.py", "patient_figures.py"]
if __name__ == "__main__":
    for s in STEPS:
        print(f"\n=== {s}")
        r = subprocess.run([sys.executable, s], cwd="src")
        if r.returncode:
            sys.exit(f"step failed: {s}")
    subprocess.run([sys.executable, "-m", "pytest", "-q", "tests"])
