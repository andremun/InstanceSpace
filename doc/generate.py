#!/usr/bin/env python3
"""
doc/generate.py -- InstanceSpace documentation site generator.

Toolchain decision recorded in llm-server/docs/instancespace_docs_site_
framework.md, section 4a (trial run 2026-09-09, confirmed by Andres):
Markdown-first source (doc/src/*.md) -> HTML (doc/html/*.html), via
python-markdown, sharing one stylesheet (doc/style.css).

Nested opts.*/out.* struct fields (the house-style convention throughout
core/ -- see ISAdefaults.m) use the DOTTED-KEY FLATTENING convention:
plain Markdown tables cannot nest a table inside a cell, so a field like
opts.dims is written as its own row (`opts.dims` in the key column,
lightly indented with a non-breaking-space prefix) rather than as a
sub-table nested inside the `opts` row. Confirmed in the trial to render
at least as consistently as a nested-table alternative across a page
with many fields.

Usage
-----
    python doc/generate.py            # regenerate every page under doc/src/
    python doc/generate.py PILOT      # regenerate just doc/src/PILOT.md
"""

import sys
from pathlib import Path

import markdown

DOC_DIR = Path(__file__).parent
SRC_DIR = DOC_DIR / "src"
OUT_DIR = DOC_DIR / "html"

NAV_DATA = [
    ("InstanceSpace", [("Landing.html", "InstanceSpace")]),
    ("Getting Started", [("GettingStarted.html", "Getting Started"), ("InteractiveWalkthrough.html", "Interactive Walkthrough")]),
    ("Standalone pipeline functions", [("PRELIM.html", "PRELIM"), ("SIFTED.html", "SIFTED"), ("PILOT.html", "PILOT"), ("CLOISTER.html", "CLOISTER"), ("PYTHIA.html", "PYTHIA"), ("TRACE.html", "TRACE"), ("FILTER.html", "FILTER"), ("PILOTviewpoint.html", "PILOTviewpoint"), ("INIT.html", "INIT")]),
    ("Class interface", [("InstanceSpace.html", "InstanceSpace")]),
    ("Backward-compatibility wrappers", [("buildIS.html", "buildIS"), ("exploreIS.html", "exploreIS")]),
    ("Utilities", [("ISAdefaults.html", "ISAdefaults"), ("ISAvalidateOpts.html", "ISAvalidateOpts"), ("ISAsubsetData.html", "ISAsubsetData"), ("ISAgetClassifierFcn.html", "ISAgetClassifierFcn"), ("ISAmigrateModel.html", "ISAmigrateModel")]),
    ("Output", [("scriptcsv.html", "scriptcsv"), ("scriptpng.html", "scriptpng"), ("scriptweb.html", "scriptweb"), ("ISArecallView.html", "ISArecallView")]),
    ("Deprecated Functions", [("Deprecated.html", "Deprecated")]),
    ("Reference", [("OptionsReference.html", "Options Reference"), ("MetadataFormat.html", "Metadata File Format"), ("MigratingLegacyModel.html", "Migrating a Legacy Model")]),
    ("What's New", [("WhatsNew.html", "What's New")]),
]

TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<link rel="stylesheet" href="../style.css">
</head>
<body>
<div class="sidebar">
{nav}
</div>
<div class="content">
{body}
<hr class="footer-rule">
<p class="copyright"><em>Copyright (c) 2026 Mario Andres Munoz Acosta and contributors,
School of Computing and Information Systems, The University of Melbourne.
Released under the PolyForm Noncommercial License 1.0.0.</em></p>
</div>
</body>
</html>
"""


def render_nav(current_stem: str) -> str:
    """Build sidebar HTML from NAV_DATA, marking current_stem's entry as <span class="current-page">."""
    nav_items = []
    for section_label, links in NAV_DATA:
        nav_items.append(f'<h4>{section_label}</h4>')
        nav_items.append('<ul>')
        for href, text in links:
            if href == current_stem + ".html":
                nav_items.append(f'<li><span class="current-page">{text}</span></li>')
            else:
                nav_items.append(f'<li><a href="{href}">{text}</a></li>')
        nav_items.append('</ul>')
    return '\n'.join(nav_items)


def render(src_path: Path) -> Path:
    """Render one doc/src/<Name>.md -> doc/html/<Name>.html."""
    text = src_path.read_text()
    lines = text.splitlines()
    if not lines or not lines[0].startswith("# "):
        raise ValueError(f"{src_path}: expected an H1 title on line 1")

    title = lines[0].lstrip("# ").strip()
    # Line 2 is blank, line 3 (if present and not itself a heading) is the
    # tagline -- matches this project's house-style docstring shape (H1
    # comment line, then a one-line description before Inputs/Outputs).
    tagline = ""
    body_start = 1
    if len(lines) > 2 and lines[1].strip() == "" and not lines[2].startswith("#"):
        tagline = lines[2].strip()
        body_start = 3

    rest = "\n".join(lines[body_start:])
    html_body = markdown.markdown(rest, extensions=["tables", "fenced_code"])
    html_body = html_body.replace("<pre><code>", '<pre class="syntax"><code>', 1)
    html_body = html_body.replace("<table>", '<table class="args">')

    tagline_html = f'<p class="tagline">{tagline}</p>\n' if tagline else ""
    full_body = f"<h1>{title}</h1>\n{tagline_html}{html_body}"

    out_path = OUT_DIR / f"{src_path.stem}.html"
    nav_html = render_nav(src_path.stem)
    full_html = TEMPLATE.format(title=title, body=full_body, nav=nav_html)
    out_path.write_text(full_html)
    return out_path


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    targets = sys.argv[1:]
    if targets:
        src_paths = [SRC_DIR / f"{name}.md" for name in targets]
    else:
        src_paths = sorted(SRC_DIR.glob("*.md"))

    if not src_paths:
        print(f"No .md sources found under {SRC_DIR}")
        return

    for src_path in src_paths:
        if not src_path.exists():
            print(f"  [SKIP] {src_path} not found")
            continue
        out_path = render(src_path)
        print(f"  {src_path.relative_to(DOC_DIR)} -> {out_path.relative_to(DOC_DIR)}")


if __name__ == "__main__":
    main()
