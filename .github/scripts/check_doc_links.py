#!/usr/bin/env python3
"""Fail if a page in doc/html links to a missing page, asset or anchor,
or repeats an element id.

Usage: python3 .github/scripts/check_doc_links.py [doc/html]
"""
import html
import re
import sys
from pathlib import Path


def main(root):
    root = Path(root)
    pages = {p.name: p.read_text(encoding="utf-8") for p in root.glob("*.html")}
    ids = {name: re.findall(r'id="([^"]+)"', text) for name, text in pages.items()}
    problems = []
    for name, text in pages.items():
        dups = sorted({i for i in ids[name] if ids[name].count(i) > 1})
        if dups:
            problems.append(f"{name}: duplicate ids {dups}")
        for href in re.findall(r'href="([^"]+)"', text):
            if href.startswith(("http://", "https://", "mailto:")):
                continue
            target, _, anchor = href.partition("#")
            target = target or name
            if target not in pages:
                if not (root / target).is_file():
                    problems.append(f"{name}: broken link {href}")
                continue
            if anchor and html.unescape(anchor) not in ids[target]:
                problems.append(f"{name}: missing anchor {href}")
    for p in problems:
        print(p)
    print(f"Checked {len(pages)} pages: {len(problems)} problem(s).")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "doc/html"))
