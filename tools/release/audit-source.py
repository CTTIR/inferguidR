#!/usr/bin/env python3
"""Inspect a source archive without extracting it or executing its contents."""

import argparse
import datetime
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        parser.error("output must be a new or empty directory")
    args.output.mkdir(parents=True, exist_ok=True)
    failures, inventory = [], []
    archive_hash = None
    description = None
    try:
        archive_hash = hashlib.sha256(args.archive.read_bytes()).hexdigest()
        with tarfile.open(args.archive, "r:gz") as source:
            members = source.getmembers()
            seen = set()
            regular_files = {}
            for member in members:
                name = member.name
                path = PurePosixPath(name)
                canonical = str(path)
                if canonical in seen:
                    failures.append(f"Duplicate archive member: {name}")
                seen.add(canonical)
                if name.rstrip("/") != canonical:
                    failures.append(f"Noncanonical archive path: {name}")
                if path.is_absolute() or ".." in path.parts or "\\" in name:
                    failures.append(f"Unsafe archive path: {name}")
                    continue
                if not path.parts or path.parts[0] != "inferguidR":
                    failures.append(f"Unexpected package root: {name}")
                    continue
                if not (member.isfile() or member.isdir()):
                    failures.append(f"Link or special archive member: {name}")
                    continue
                relative = path.parts[1:]
                if relative and relative[0] in {
                    ".git", ".github", ".local-private", "project", "specs", "validation",
                    "renv", "admin", "AGENTS.md", "CLAUDE.md", "ROADMAP.md", "VALIDATION.md",
                }:
                    failures.append(f"Development/private material in package: {name}")
                if path.suffix.lower() in {".exe", ".dll", ".so", ".dylib", ".ttf", ".otf"}:
                    failures.append(f"Unexpected binary or font in pure-R source package: {name}")
                if not member.isfile():
                    continue
                if member.size > 32 * 1024 * 1024:
                    failures.append(f"Oversized member requires explicit review: {name}")
                    continue
                payload = source.extractfile(member).read()
                regular_files[canonical] = len(payload)
                inventory.append({"path": name, "bytes": len(payload),
                                  "sha256": hashlib.sha256(payload).hexdigest()})
                if name == "inferguidR/DESCRIPTION":
                    description = payload.decode("utf-8")
                if path.suffix.lower() in {".r", ".rd", ".rmd", ".md", ".json", ".html", ".txt", ".yml", ".yaml"}:
                    for pattern in (rb"/(?:home|Users)/[^/\s]+/", rb"/data/GitHome/",
                                    rb"/data/GitHub/[^/\s]+/private/"):
                        if re.search(pattern, payload):
                            failures.append(f"Private/runtime path reference in {name}")
            for required in ("DESCRIPTION", "NAMESPACE", "LICENSE"):
                if not regular_files.get(f"inferguidR/{required}", 0):
                    failures.append(f"Required package file is absent, empty or not regular: {required}")
        if description is not None:
            if not re.search(r"^Package:\s*inferguidR\s*$", description, re.M):
                failures.append("DESCRIPTION package identity mismatch")
            if not re.search(r"^License:\s*MIT \+ file LICENSE\s*$", description, re.M):
                failures.append("DESCRIPTION licence requires review")
    except (OSError, tarfile.TarError, UnicodeError) as error:
        failures.append(str(error))
    receipt = {
        "scope": "source archive identity, inventory and distribution hygiene; not full licence or release certification",
        "checked_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "archive": str(args.archive.resolve()), "archive_sha256": archive_hash,
        "file_count": len(inventory), "uncompressed_file_bytes": sum(x["bytes"] for x in inventory),
        "status": "FAIL" if failures else "PASS", "failures": failures,
    }
    (args.output / "inventory.json").write_text(json.dumps(inventory, indent=2) + "\n")
    (args.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(f"Source archive audit: {receipt['status']}; {len(inventory)} files; {len(failures)} findings")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
