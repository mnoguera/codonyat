from pathlib import Path
import sys

# Ensure tests import the local package even without installing it into the venv
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))