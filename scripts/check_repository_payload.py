#!/usr/bin/env python3
"""Fail if a Git commit contains local run artifacts or oversized files."""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAX_BYTES = 10 * 1024 * 1024
FORBIDDEN_ROOTS = {"release", "runs", "work", "inputs", "analysis", "tools"}
FORBIDDEN_SUFFIXES = {
    ".a3m", ".cif", ".pdb", ".npz", ".npy", ".parquet", ".pt", ".pth",
    ".ckpt", ".safetensors", ".log", ".lock", ".metadata", ".pyc",
}


def git_paths(mode: str) -> list[str]:
    if mode == "staged":
        command = ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z"]
    else:
        command = ["git", "ls-files", "-z"]
    output = subprocess.check_output(command, cwd=ROOT)
    return [item.decode() for item in output.split(b"\0") if item]


def forbidden_reason(relative: str, path: Path) -> str | None:
    parts = Path(relative).parts
    if not parts:
        return None
    if parts[0] in FORBIDDEN_ROOTS:
        return f"local/generated root: {parts[0]}/"
    if relative.startswith("intro/.bayesian_preregistration_history/"):
        return "editor history"
    if path.is_symlink():
        target = path.readlink()
        if target.is_absolute() or ".." in target.parts:
            return f"external symlink: {target}"
    lower = relative.lower()
    if lower.endswith((".cif.gz", ".pdb.gz")):
        return "coordinate archive"
    if path.suffix.lower() in FORBIDDEN_SUFFIXES:
        return f"forbidden artifact extension: {path.suffix.lower()}"
    if path.is_file() and path.stat().st_size > MAX_BYTES:
        return f"file exceeds {MAX_BYTES} bytes"
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("staged", "all"), default="staged")
    args = parser.parse_args()
    failures: list[tuple[str, str]] = []
    for relative in git_paths(args.mode):
        path = ROOT / relative
        reason = forbidden_reason(relative, path)
        if reason:
            failures.append((relative, reason))
    if failures:
        print("Repository payload check: FAIL")
        for relative, reason in failures:
            print(f"- {relative}: {reason}")
        return 1
    print(f"Repository payload check: PASS ({len(git_paths(args.mode))} files checked)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
