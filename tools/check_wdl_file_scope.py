#!/usr/bin/env python3
"""Reject workflow-scope file creation that needs a local Cromwell filesystem.

Parse local sources without resolving imports. Pass the full workflows directory
to cover all repository imports. Task declarations and commands are permitted.
Requires MiniWDL, already installed by the workflow validation CI job.
"""
import argparse
import sys
from pathlib import Path

import WDL


def workflow_file_writes(node):
    """Walk workflow syntax, including call inputs but not the called task."""
    if node is None:
        return []
    writes = []
    if isinstance(node, WDL.Expr.Apply) and node.function_name.startswith("write_"):
        writes.append(f"{node.function_name} at line {node.pos.line}")
    for child in node.children:
        writes.extend(workflow_file_writes(child))
    return writes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="+", type=Path, help="WDL files or directories")
    arguments = parser.parse_args()
    files = set()
    for path in arguments.paths:
        if not path.exists():
            parser.error(f"Path does not exist: {path}")
        files.update(path.rglob("*.wdl") if path.is_dir() else [path])
    if not files:
        parser.error("No WDL files found")
    failures = []
    for path in sorted(files):
        try:
            document = WDL.parse_document(path.read_text(), uri=str(path))
            failures.extend(f"{path}: {issue}" for issue in workflow_file_writes(document.workflow))
        except (OSError, WDL.Error.SyntaxError) as error:
            failures.append(f"{path}: cannot parse WDL: {error}")
    if failures:
        print("Terra file-scope check failed. Create files in task scope:", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1
    print(f"Terra file-scope check passed for {len(files)} WDL files.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
