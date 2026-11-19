"""Reassemble and verify the exact R2 reproducibility ZIP (Python 3, standard library)."""
import argparse
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parent
meta = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=Path(meta["name"]))
out = parser.parse_args().output
total_hash = hashlib.sha256()
size = 0
with out.open("xb") as f:
    try:
        for part in meta["parts"]:
            data = (root / "parts" / part["name"]).read_bytes()
            if len(data) != part["size_bytes"] or hashlib.sha256(data).hexdigest() != part["sha256"]:
                raise ValueError(f"Part failed verification: {part['name']}")
            f.write(data)
            total_hash.update(data)
            size += len(data)
        if size != meta["size_bytes"] or total_hash.hexdigest() != meta["sha256"]:
            raise ValueError("Final ZIP failed verification")
    except BaseException:
        f.close()
        out.unlink()
        raise
print(f"Verified: {out} ({size:,} bytes)")
