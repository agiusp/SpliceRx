import sys
from pathlib import Path

# app/backend/ — so `import sjlookup` resolves to the combined-app package.
BACKEND = Path(__file__).resolve().parents[2]
if str(BACKEND) not in sys.path:
    sys.path.insert(0, str(BACKEND))
