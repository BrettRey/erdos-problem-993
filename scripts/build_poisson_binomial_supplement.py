#!/usr/bin/env python3
"""Build the deterministic Poisson--binomial certificate supplement."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile
from typing import Final
from zipfile import ZIP_STORED, ZipFile, ZipInfo


REPOSITORY_ROOT: Final = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT: Final = (
    REPOSITORY_ROOT
    / "paper"
    / "poisson_binomial"
    / "poisson_binomial_certificate_supplement.zip"
)
FIXED_TIMESTAMP: Final = (1980, 1, 1, 0, 0, 0)
REGULAR_FILE_MODE: Final = 0o100644 << 16
LEAN_PROJECT_ROOT: Final = (
    REPOSITORY_ROOT / "formalization" / "pb_effective_drop_aristotle"
)

ARCHIVE_MEMBERS: Final = (
    (
        "CERTIFICATE.md",
        REPOSITORY_ROOT / "paper" / "poisson_binomial" / "CERTIFICATE.md",
    ),
    ("LICENSE", REPOSITORY_ROOT / "LICENSE"),
    (
        "formalization/pb_effective_drop_aristotle/ARISTOTLE_RESULT.md",
        LEAN_PROJECT_ROOT / "ARISTOTLE_RESULT.md",
    ),
    (
        "formalization/pb_effective_drop_aristotle/PBReserve.lean",
        LEAN_PROJECT_ROOT / "PBReserve.lean",
    ),
    (
        "formalization/pb_effective_drop_aristotle/PBReserve/Core.lean",
        LEAN_PROJECT_ROOT / "PBReserve" / "Core.lean",
    ),
    (
        "formalization/pb_effective_drop_aristotle/PROMPT.md",
        LEAN_PROJECT_ROOT / "PROMPT.md",
    ),
    (
        "formalization/pb_effective_drop_aristotle/PROOF_CONTEXT.md",
        LEAN_PROJECT_ROOT / "PROOF_CONTEXT.md",
    ),
    (
        "formalization/pb_effective_drop_aristotle/lake-manifest.json",
        LEAN_PROJECT_ROOT / "lake-manifest.json",
    ),
    (
        "formalization/pb_effective_drop_aristotle/lakefile.toml",
        LEAN_PROJECT_ROOT / "lakefile.toml",
    ),
    (
        "formalization/pb_effective_drop_aristotle/lean-toolchain",
        LEAN_PROJECT_ROOT / "lean-toolchain",
    ),
    (
        "results/universal_pb_finite_bernstein_certificate_2026-07-10.json",
        REPOSITORY_ROOT
        / "results"
        / "universal_pb_finite_bernstein_certificate_2026-07-10.json",
    ),
    (
        "results/universal_pb_finite_bernstein_full_certificate_2026-07-16.json",
        REPOSITORY_ROOT
        / "results"
        / "universal_pb_finite_bernstein_full_certificate_2026-07-16.json",
    ),
    (
        "scripts/build_poisson_binomial_supplement.py",
        REPOSITORY_ROOT / "scripts" / "build_poisson_binomial_supplement.py",
    ),
    (
        "scripts/check_universal_pb_finite_bernstein_certificate.py",
        REPOSITORY_ROOT
        / "scripts"
        / "check_universal_pb_finite_bernstein_certificate.py",
    ),
    (
        "scripts/verify_universal_pb_finite_bernstein.py",
        REPOSITORY_ROOT
        / "scripts"
        / "verify_universal_pb_finite_bernstein.py",
    ),
)


def digest(data: bytes) -> str:
    """Return the hexadecimal SHA-256 digest of *data*."""
    return sha256(data).hexdigest()


def zip_info(name: str) -> ZipInfo:
    """Return normalized metadata for a deterministic regular-file member."""
    info = ZipInfo(name, FIXED_TIMESTAMP)
    info.compress_type = ZIP_STORED
    info.create_system = 3
    info.external_attr = REGULAR_FILE_MODE
    return info


def source_payloads() -> dict[str, bytes]:
    """Read and validate all source members in canonical archive order."""
    payloads: dict[str, bytes] = {}
    for name, path in sorted(ARCHIVE_MEMBERS):
        if not path.is_file():
            raise FileNotFoundError(f"required supplement member is missing: {path}")
        payloads[name] = path.read_bytes()
    return payloads


def manifest_bytes(payloads: dict[str, bytes]) -> bytes:
    """Return a standard SHA-256 manifest for all non-manifest members."""
    lines = [f"{digest(data)}  {name}\n" for name, data in payloads.items()]
    return "".join(lines).encode("ascii")


def verify_archive(path: Path, expected: dict[str, bytes]) -> None:
    """Reopen *path* and require exact member names, metadata, and bytes."""
    with ZipFile(path, "r") as archive:
        names = archive.namelist()
        if names != list(expected):
            raise RuntimeError(f"unexpected ZIP member order: {names}")
        if archive.testzip() is not None:
            raise RuntimeError("ZIP integrity check failed")
        for name, data in expected.items():
            info = archive.getinfo(name)
            if info.date_time != FIXED_TIMESTAMP:
                raise RuntimeError(f"noncanonical timestamp for {name}")
            if info.compress_type != ZIP_STORED:
                raise RuntimeError(f"noncanonical compression for {name}")
            if info.external_attr != REGULAR_FILE_MODE:
                raise RuntimeError(f"noncanonical file mode for {name}")
            if archive.read(name) != data:
                raise RuntimeError(f"content mismatch for {name}")


def build_archive(output: Path) -> dict[str, object]:
    """Build, verify, and checksum the supplement at *output*."""
    payloads = source_payloads()
    payloads["MANIFEST.sha256"] = manifest_bytes(payloads)
    payloads = dict(sorted(payloads.items()))

    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        with ZipFile(temporary_path, "w", compression=ZIP_STORED) as archive:
            for name, data in payloads.items():
                archive.writestr(zip_info(name), data)
        verify_archive(temporary_path, payloads)
        os.replace(temporary_path, output)
        output.chmod(0o644)
        temporary_path = None
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)

    archive_digest = digest(output.read_bytes())
    checksum_path = output.with_suffix(output.suffix + ".sha256")
    checksum_path.write_text(
        f"{archive_digest}  {output.name}\n",
        encoding="ascii",
    )
    return {
        "status": "passed",
        "archive": str(output),
        "archive_sha256": archive_digest,
        "checksum": str(checksum_path),
        "members": len(payloads),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Destination ZIP path.",
    )
    args = parser.parse_args()
    print(json.dumps(build_archive(args.out), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
