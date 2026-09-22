"""Repository paths and the relocation map for historical source manifests."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCATIONS = json.loads((ROOT / 'configs/file_locations.json').read_text())


def project_file(name):
    """Resolve an old root-level filename without changing historical manifests."""
    return ROOT / LOCATIONS.get(str(name), str(name))
