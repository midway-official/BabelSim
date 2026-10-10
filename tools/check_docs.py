#!/usr/bin/env python3
"""Check maintained Markdown links, heading anchors, and documented Make targets.

Frozen experiment source copies are inputs, not maintained manuals. Historical
links to unavailable development artifacts are reported separately; the archive
index and all current documentation must resolve without exceptions.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import unicodedata
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"!?\[[^\]\n]*\]\(\s*(<[^>\n]+>|[^\s)]+)(?:\s+[\"'][^\n]*?[\"'])?\s*\)")


def headings(path: Path) -> set[str]:
    result: set[str] = set()
    counts: dict[str, int] = {}
    fenced = False
    for line in path.read_text().splitlines():
        if re.match(r"^\s*(```|~~~)", line):
            fenced = not fenced
            continue
        match = re.match(r"^#{1,6}\s+(.+?)\s*#*\s*$", line)
        if fenced or not match:
            continue
        title = re.sub(r"<[^>]+>", "", match.group(1)).lower()
        slug = "".join(c for c in title if c in " -_" or unicodedata.category(c)[0] in "LN")
        slug = slug.replace(" ", "-")
        count = counts.get(slug, 0)
        counts[slug] = count + 1
        result.add(slug if count == 0 else f"{slug}-{count}")
    return result


def documents(delivery_only: bool = False) -> list[Path]:
    selected = [ROOT / "README.md", ROOT / "dependencies/README.md"]
    selected.extend((ROOT / "docs").glob("*.md"))
    for group in ("user-guide", "reference", "developer-guide", "examples"):
        selected.extend((ROOT / "docs" / group).rglob("*.md"))
    if not delivery_only:
        selected.extend((ROOT / "archive/docs").glob("*.md"))
        selected.extend((ROOT / "archive/docs/manuals").glob("*.md"))
        selected.extend((ROOT / "archive/docs/reports").glob("*.md"))
        selected.extend((ROOT / "archive/docs/reports").glob("*-evidence/README.md"))
        selected.extend((ROOT / "archive/docs/history").rglob("*.md"))
    # Only tracked case entry documents; leave unrelated local cases alone.
    if not delivery_only:
        names = subprocess.check_output(
            ["git", "ls-files", "-z", "cases/**/README.md"], cwd=ROOT
        ).decode().split("\0")
        selected.extend(ROOT / name for name in names if name)
    return sorted(set(p for p in selected if p.is_file()))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, help="optional diagnostic output")
    parser.add_argument("--delivery-only", action="store_true", help="check product manuals without Git or optional archives")
    args = parser.parse_args()
    makefile = (ROOT / "Makefile").read_text()
    targets = set(re.findall(r"^([a-z][a-z0-9-]*):", makefile, re.MULTILINE))
    errors: list[dict[str, object]] = []
    historical: list[dict[str, object]] = []
    checked_links = 0
    anchors: dict[Path, set[str]] = {}
    navigation: dict[Path, set[Path]] = {}
    docs = documents(args.delivery_only)
    for path in docs:
        relative = path.relative_to(ROOT).as_posix()
        is_historical = relative.startswith(("archive/docs/history/reports/", "archive/docs/history/plans/"))
        content = path.read_text()
        # Code fences contain examples, not rendered Markdown links.
        rendered = re.sub(r"(?ms)^\s*(```|~~~).*?^\s*\1\s*$", "", content)
        for match in LINK.finditer(rendered):
            destination = match.group(1).strip("<>")
            parsed = urlsplit(destination)
            if parsed.scheme or parsed.netloc:
                continue
            checked_links += 1
            target = (path.parent / unquote(parsed.path)).resolve() if parsed.path else path
            if target.suffix == ".md":
                navigation.setdefault(path, set()).add(target)
            problem = None
            if not target.exists():
                problem = "missing file or directory"
            elif parsed.fragment and target.suffix == ".md":
                if target not in anchors:
                    anchors[target] = headings(target)
                if unquote(parsed.fragment) not in anchors[target]:
                    problem = "missing heading anchor"
            if problem:
                item = dict(file=relative, target=destination, problem=problem)
                (historical if is_historical else errors).append(item)
        if not is_historical:
            for target in set(re.findall(r"\bmake\s+(?:-[^\s`]+\s+)*((?:test|validate)(?:-[a-z0-9]+)*)\b", content)):
                if target not in targets:
                    errors.append(dict(file=relative, target=target, problem="missing Make target"))
        if is_historical and "历史开发记录" not in content:
            errors.append(dict(file=relative, problem="missing historical notice"))
    # Every delivery manual must be discoverable from the product docs entrance.
    entrance = ROOT / "docs/README.md"
    reached: set[Path] = set()
    pending = [entrance]
    while pending:
        path = pending.pop()
        if path in reached:
            continue
        reached.add(path)
        pending.extend(target for target in navigation.get(path, set()) - reached
                       if target.is_relative_to(ROOT / "docs"))
    manuals = [path for path in docs if path.is_relative_to(ROOT / "docs")]
    for path in manuals:
        if path not in reached:
            errors.append(dict(file=path.relative_to(ROOT).as_posix(), problem="manual is unreachable from docs/README.md"))
    allowed_groups = {"user-guide", "reference", "developer-guide", "examples"}
    for path in (ROOT / "docs").iterdir():
        if (path.is_dir() and path.name not in allowed_groups) or (path.is_file() and path.name != "README.md"):
            errors.append(dict(file=path.relative_to(ROOT).as_posix(), problem="unexpected delivery docs entry"))
    result = dict(status="passed" if not errors else "failed", documents=len(docs),
                  localLinks=checked_links, errors=errors, historicalUnresolved=historical,
                  deliveryManuals=len(manuals), deliveryManualsReachable=sum(path in reached for path in manuals),
                  deliveryOnly=args.delivery_only,
                  frozenSourceSnapshotsExcluded=True)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
