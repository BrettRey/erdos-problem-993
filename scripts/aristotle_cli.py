#!/usr/bin/env python3
"""Thin repo-local wrapper around the official Aristotle CLI."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DOWNLOADS = Path.home() / "Downloads"


def find_aristotle_bin() -> str:
    override = os.environ.get("ARISTOTLE_BIN")
    if override:
        return override
    found = shutil.which("aristotle")
    if found:
        return found
    local = Path.home() / ".local" / "bin" / "aristotle"
    if local.exists():
        return str(local)
    raise SystemExit(
        "Aristotle CLI not found. Install it with:\n"
        "  uv tool install aristotlelib"
    )


def require_api_key() -> None:
    if os.environ.get("ARISTOTLE_API_KEY"):
        return
    raise SystemExit(
        "ARISTOTLE_API_KEY is not set.\n"
        "Create one in Aristotle Dashboard -> API Keys, then export it, e.g.:\n"
        '  echo \'export ARISTOTLE_API_KEY=\"...\"\' >> ~/.zshrc\n'
        "  source ~/.zshrc"
    )


def timestamped_destination(prefix: str) -> Path:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return DEFAULT_DOWNLOADS / f"{prefix}-{stamp}.tar.gz"


def run(cmd: list[str]) -> int:
    completed = subprocess.run(cmd, cwd=REPO_ROOT)
    return completed.returncode


def resolve_prompt(args: argparse.Namespace) -> str:
    if args.prompt_file:
        content = args.prompt_file.read_text(encoding="utf-8").strip()
        if not content:
            raise SystemExit(f"Prompt file is empty: {args.prompt_file}")
        return content
    return args.prompt


def cmd_submit(args: argparse.Namespace) -> int:
    require_api_key()
    prompt = resolve_prompt(args)
    cmd = [
        find_aristotle_bin(),
        "submit",
        prompt,
        "--project-dir",
        str(args.project_dir.resolve()),
    ]
    if args.wait:
        cmd.append("--wait")
    if args.destination is not None:
        cmd.extend(["--destination", str(args.destination.resolve())])
    elif args.wait:
        cmd.extend(
            ["--destination", str(timestamped_destination("aristotle-submit"))]
        )
    return run(cmd)


def cmd_formalize(args: argparse.Namespace) -> int:
    require_api_key()
    cmd = [find_aristotle_bin(), "formalize", str(args.input_file.resolve())]
    if args.wait:
        cmd.append("--wait")
    if args.destination is not None:
        cmd.extend(["--destination", str(args.destination.resolve())])
    elif args.wait:
        default_name = f"{args.input_file.stem}-aristotle"
        cmd.extend(["--destination", str(timestamped_destination(default_name))])
    return run(cmd)


def cmd_result(args: argparse.Namespace) -> int:
    require_api_key()
    cmd = [find_aristotle_bin(), "result", args.project_id]
    if args.wait:
        cmd.append("--wait")
    if args.destination is not None:
        cmd.extend(["--destination", str(args.destination.resolve())])
    elif args.wait:
        cmd.extend(
            ["--destination", str(DEFAULT_DOWNLOADS / f"{args.project_id}-aristotle.tar.gz")]
        )
    return run(cmd)


def cmd_list(args: argparse.Namespace) -> int:
    require_api_key()
    cmd = [find_aristotle_bin(), "list", "--limit", str(args.limit)]
    if args.pagination_key:
        cmd.extend(["--pagination-key", args.pagination_key])
    for status in args.status:
        cmd.extend(["--status", status])
    return run(cmd)


def cmd_cancel(args: argparse.Namespace) -> int:
    require_api_key()
    return run([find_aristotle_bin(), "cancel", args.project_id])


def cmd_doctor(_: argparse.Namespace) -> int:
    bin_path = find_aristotle_bin()
    print(f"repo_root={REPO_ROOT}")
    print(f"aristotle_bin={bin_path}")
    print(f"api_key={'set' if os.environ.get('ARISTOTLE_API_KEY') else 'missing'}")
    print(f"downloads_dir={DEFAULT_DOWNLOADS}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Repo-local wrapper around the official Aristotle CLI."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    submit_parser = subparsers.add_parser(
        "submit", help="Submit the current Lean project or another project dir."
    )
    submit_group = submit_parser.add_mutually_exclusive_group(required=True)
    submit_group.add_argument("prompt", nargs="?", help="Prompt text for Aristotle.")
    submit_group.add_argument(
        "--prompt-file",
        type=Path,
        help="Read the prompt from a file, e.g. aristotle_prompt_stp2_closure.md.",
    )
    submit_parser.add_argument(
        "--project-dir",
        type=Path,
        default=REPO_ROOT,
        help=f"Project directory to upload. Default: {REPO_ROOT}",
    )
    submit_parser.add_argument(
        "--wait", action="store_true", help="Wait for completion and download the result."
    )
    submit_parser.add_argument(
        "--destination",
        type=Path,
        help="Output tarball path. If omitted with --wait, a timestamped file goes to ~/Downloads.",
    )
    submit_parser.set_defaults(func=cmd_submit)

    formalize_parser = subparsers.add_parser(
        "formalize", help="Formalize a single input file with Aristotle."
    )
    formalize_parser.add_argument("input_file", type=Path, help="File to formalize.")
    formalize_parser.add_argument("--wait", action="store_true")
    formalize_parser.add_argument("--destination", type=Path)
    formalize_parser.set_defaults(func=cmd_formalize)

    result_parser = subparsers.add_parser(
        "result", help="Fetch a result for an existing Aristotle project id."
    )
    result_parser.add_argument("project_id")
    result_parser.add_argument("--wait", action="store_true")
    result_parser.add_argument("--destination", type=Path)
    result_parser.set_defaults(func=cmd_result)

    list_parser = subparsers.add_parser("list", help="List recent Aristotle jobs.")
    list_parser.add_argument("--limit", type=int, default=10)
    list_parser.add_argument("--pagination-key")
    list_parser.add_argument(
        "--status",
        nargs="*",
        default=[],
        choices=[
            "NOT_STARTED",
            "QUEUED",
            "IN_PROGRESS",
            "COMPLETE",
            "COMPLETE_WITH_ERRORS",
            "OUT_OF_BUDGET",
            "FAILED",
            "CANCELED",
        ],
    )
    list_parser.set_defaults(func=cmd_list)

    cancel_parser = subparsers.add_parser("cancel", help="Cancel a queued or running job.")
    cancel_parser.add_argument("project_id")
    cancel_parser.set_defaults(func=cmd_cancel)

    doctor_parser = subparsers.add_parser(
        "doctor", help="Show local Aristotle CLI and auth configuration."
    )
    doctor_parser.set_defaults(func=cmd_doctor)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
