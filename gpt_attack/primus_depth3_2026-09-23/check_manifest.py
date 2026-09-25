#!/usr/bin/env python3
"""Verify every file listed in the SHA-256 manifest, without third-party imports."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath


def main() -> None:
    root = Path(__file__).resolve().parent
    seen: set[str] = set()
    for line in (root / "MANIFEST.sha256").read_text().splitlines():
        digest, name = line.split("  ", 1)
        path = PurePosixPath(name)
        if (len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest)
                or path.is_absolute() or ".." in path.parts or name in seen):
            raise ValueError(f"Invalid or duplicate manifest entry: {name}")
        target = root.joinpath(*path.parts)
        if not target.resolve().is_relative_to(root) or not target.is_file():
            raise ValueError(f"Missing or escaping file: {name}")
        actual = hashlib.sha256(target.read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f"SHA-256 mismatch: {name}")
        seen.add(name)
    if not seen:
        raise ValueError("Empty manifest")
    print(json.dumps({"status": "PASS", "files_verified": len(seen)}, indent=2))


if __name__ == "__main__":
    main()
