"""Make every project folder importable so tests can `import three_body`, etc."""
import sys
from pathlib import Path

for project in sorted((Path(__file__).resolve().parents[1] / "projects").iterdir()):
    if project.is_dir():
        sys.path.insert(0, str(project))
