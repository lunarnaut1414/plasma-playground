"""Make the repo root importable so `import plasmaplay` works in tests
without requiring an editable install. (pytest imports this automatically.)"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
