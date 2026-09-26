#!/usr/bin/env python3
"""Build the InstanceSpace reference documentation from doc/src/*.md.

Writes doc/html/, which serves two readers from the same files:
  * MATLAB's Help browser, through info.xml (repo root) and the
    doc/html/helptoc.xml written here.
  * GitHub Pages, through .github/workflows/docs-pages.yml.

The page layout follows the MathWorks reference-page convention: a purpose
line under the title, a Syntax block linked to its Description paragraphs,
collapsible Examples, expandable Input/Output Arguments, Version History,
and See Also. The styling is an original look-alike (doc/assets/); no
MathWorks CSS, fonts, or images are used.

Source conventions (doc/src/<Page>.md)
-------------------------------------
  # Title                     line 1
  One-line purpose            line 3 (after one blank line)
  <!-- opts: pilot -->        optional: resolves bare `opts.x` code spans
                              on this page to `opts.pilot.x`
  ## Syntax                   fenced block, one calling form per line
  ## Description              a paragraph that opens with a calling form in
                              backticks is linked from the Syntax block
  ## Examples                 each ### is one collapsible example
  ## Input Arguments          each ### is one expandable argument, written
  ## Output Arguments         as  ### `name` — Short description
  ## Properties               (#### nests fields inside an argument)
  ## Version History          each ### is one release entry

Code spans that name a documented page (`PILOT`) or an option
(`opts.pilot.alpha`) become links automatically.

Usage
-----
    python doc/generate.py           # rebuild everything
    python doc/generate.py --check   # fail if doc/html is out of date
"""

import html
import json
import re
import sys
from pathlib import Path

import markdown

DOC_DIR = Path(__file__).resolve().parent
REPO_DIR = DOC_DIR.parent
SRC_DIR = DOC_DIR / "src"
ASSET_DIR = DOC_DIR / "assets"
OUT_DIR = DOC_DIR / "html"

TOOLBOX = "Instance Space Analysis Toolbox"
REPO_URL = "https://github.com/andremun/InstanceSpace"

# Table of contents: (section label, landing page, [(page stem, label)]).
# The single source for the sidebar, breadcrumbs, the Functions page and
# helptoc.xml, so the Help browser and the web site cannot drift apart.
TOC = [
    ("Getting Started", "GettingStarted", [
        ("GettingStarted", "Getting Started with Instance Space Analysis"),
        ("InteractiveWalkthrough", "Stage-by-Stage Walkthrough"),
        ("MetadataFormat", "Metadata File Format"),
    ]),
    ("Pipeline Stages", "FunctionList", [
        ("INIT", "INIT"),
        ("PRELIM", "PRELIM"),
        ("FILTER", "FILTER"),
        ("SIFTED", "SIFTED"),
        ("PILOT", "PILOT"),
        ("PILOTviewpoint", "PILOTviewpoint"),
        ("CLOISTER", "CLOISTER"),
        ("PYTHIA", "PYTHIA"),
        ("TRACE", "TRACE"),
    ]),
    ("Class Interface", "InstanceSpace", [
        ("InstanceSpace", "InstanceSpace"),
    ]),
    ("Backward-Compatible Wrappers", "buildIS", [
        ("buildIS", "buildIS"),
        ("exploreIS", "exploreIS"),
    ]),
    ("Utilities", "FunctionList", [
        ("ISAdefaults", "ISAdefaults"),
        ("ISAvalidateOpts", "ISAvalidateOpts"),
        ("ISAgetClassifierFcn", "ISAgetClassifierFcn"),
        ("ISAmigrateModel", "ISAmigrateModel"),
        ("ISAsubsetData", "ISAsubsetData"),
    ]),
    ("Output", "FunctionList", [
        ("scriptcsv", "scriptcsv"),
        ("scriptpng", "scriptpng"),
        ("scriptweb", "scriptweb"),
        ("ISArecallView", "ISArecallView"),
    ]),
    ("Topics", "OptionsReference", [
        ("OptionsReference", "Options Reference"),
        ("MigratingLegacyModel", "Migrating a Legacy Model"),
        ("Deprecated", "Deprecated Functions"),
    ]),
    ("Release Notes", "WhatsNew", [
        ("WhatsNew", "What's New"),
    ]),
]

# Sections whose ### headings become expandable argument blocks.
ARG_SECTIONS = {"Input Arguments", "Output Arguments", "Properties",
                "Name-Value Arguments"}

MATLAB_KEYWORDS = {
    "break", "case", "catch", "classdef", "continue", "else", "elseif",
    "end", "for", "function", "global", "if", "methods", "otherwise",
    "parfor", "persistent", "properties", "return", "switch", "try",
    "while", "true", "false",
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def slug(text):
    """Lower-case anchor id from heading or code text."""
    text = re.sub(r"<[^>]+>", "", text)
    text = html.unescape(text)
    text = re.sub(r"[^A-Za-z0-9]+", "-", text).strip("-").lower()
    return text or "section"


def strip_tags(text):
    return html.unescape(re.sub(r"<[^>]+>", " ", text))


def highlight_matlab(code):
    """Token-level MATLAB highlighting to HTML spans (no JS required)."""
    out = []
    token = re.compile(
        r"(?P<com>%.*$)"
        r"|(?P<str>\"[^\"\n]*\"|(?<![\w\)\]\}\.'])'[^'\n]*')"
        r"|(?P<num>\b\d+(?:\.\d*)?(?:[eE][+-]?\d+)?\b)"
        r"|(?P<word>\b[A-Za-z_]\w*\b)"
        r"|(?P<other>[\s\S])",
        re.M)
    for m in token.finditer(code):
        kind = m.lastgroup
        text = html.escape(m.group(0), quote=False)
        if kind == "word":
            if m.group(0) in MATLAB_KEYWORDS:
                out.append(f'<span class="kw">{text}</span>')
            else:
                out.append(text)
        elif kind == "other":
            out.append(text)
        else:
            out.append(f'<span class="{kind}">{text}</span>')
    return "".join(out)


# ---------------------------------------------------------------------------
# Page model
# ---------------------------------------------------------------------------
class Page:
    def __init__(self, path):
        self.path = path
        self.stem = path.stem
        text = path.read_text(encoding="utf-8")
        lines = text.splitlines()
        if not lines or not lines[0].startswith("# "):
            raise ValueError(f"{path}: line 1 must be an H1 title")
        self.title = lines[0][2:].strip()
        self.purpose = ""
        body_start = 1
        if len(lines) > 2 and not lines[1].strip() and lines[2].strip() \
                and not lines[2].startswith(("#", "<!--", "```", "|", "-")):
            self.purpose = lines[2].strip()
            body_start = 3
        body = "\n".join(lines[body_start:])
        meta = re.search(r"<!--\s*opts:\s*([\w.]+)\s*-->", body)
        self.opts_section = meta.group(1) if meta else None
        self.body = re.sub(r"<!--.*?-->\n?", "", body, flags=re.S)
        self.out_name = "index.html" if self.stem == "index" else f"{self.stem}.html"
        self.is_function = bool(re.search(r"^## Syntax", self.body, re.M))


def load_pages():
    pages = {}
    for path in sorted(SRC_DIR.glob("*.md")):
        page = Page(path)
        pages[page.stem] = page
    listed = {"index"} | {stem for _, _, items in TOC for stem, _ in items}
    listed.add("FunctionList")
    missing = [s for s in listed if s not in pages and s != "FunctionList"]
    unlisted = [s for s in pages if s not in listed]
    if missing or unlisted:
        raise SystemExit(f"TOC/source mismatch. Missing sources: {missing}; "
                         f"sources not in TOC: {unlisted}")
    return pages


def option_ids(pages):
    """Every documented option (opts.section and opts.section.field) -> anchor."""
    ids = {}
    ref = pages["OptionsReference"].body
    section = None
    for line in ref.splitlines():
        m = re.match(r"^## `?(opts\.\w+)`?", line)
        if m:
            section = m.group(1)
            ids[section] = slug(section)
            continue
        m = re.match(r"^\|\s*`(\w+)`\s*\|", line)
        if section and m:
            full = f"{section}.{m.group(1)}"
            ids[full] = slug(full)
    return ids


# ---------------------------------------------------------------------------
# Markdown -> HTML with the reference-page transforms
# ---------------------------------------------------------------------------
class Renderer:
    def __init__(self, pages):
        self.pages = pages
        self.opt_ids = option_ids(pages)
        self.page_names = {p.stem: p.out_name for p in pages.values()
                           if p.stem not in ("index",)}
        self.section_of = {}
        for label, _, items in TOC:
            for stem, _ in items:
                self.section_of.setdefault(stem, label)

    # -- inline links -----------------------------------------------------
    def link_code(self, html_text, page):
        """Link <code> spans naming a page or a documented option."""
        def repl(m):
            inner = m.group(1)
            raw = html.unescape(inner)
            target = None
            name = re.match(r"^([A-Za-z]\w*)(\(.*\))?$", raw)
            if name and name.group(1) in self.page_names and name.group(1) != page.stem:
                target = self.page_names[name.group(1)]
            elif raw.startswith("opts."):
                key = raw.split("=")[0].strip()
                key = re.sub(r"[^\w.].*$", "", key)
                full = key
                if full not in self.opt_ids and page.opts_section:
                    full = f"opts.{page.opts_section}.{key[5:]}"
                while full not in self.opt_ids and full.count(".") > 1:
                    full = full.rsplit(".", 1)[0]
                if full in self.opt_ids and page.stem != "OptionsReference":
                    target = f"OptionsReference.html#{self.opt_ids[full]}"
            if target:
                return f'<a class="code-link" href="{target}"><code>{inner}</code></a>'
            return m.group(0)
        # Leave code inside <pre>, headings, summaries and existing links alone.
        parts = re.split(r"(<pre.*?</pre>|<a .*?</a>|<h\d.*?</h\d>|<summary.*?</summary>)",
                         html_text, flags=re.S)
        for i in range(0, len(parts), 2):
            parts[i] = re.sub(r"<code>(.*?)</code>", repl, parts[i])
        return "".join(parts)

    # -- fenced code --------------------------------------------------------
    def md(self, text):
        """Markdown to HTML, with MATLAB-highlighted fenced code blocks."""
        blocks = []

        def stash(m):
            lang = (m.group(1) or "").strip().lower()
            code = m.group(2).rstrip("\n")
            if lang in ("matlab", "m"):
                body = highlight_matlab(code)
                cls = "code matlab"
            else:
                body = html.escape(code, quote=False)
                cls = "code"
            blocks.append(
                f'<div class="{cls}"><button class="copy" type="button" '
                f'aria-label="Copy code">Copy</button><pre><code>{body}</code></pre></div>')
            return f"\n\nCODEBLOCK{len(blocks) - 1}XX\n\n"

        text = re.sub(r"^```(\w*)[^\n]*\n(.*?)^```\s*$", stash, text, flags=re.S | re.M)
        out = markdown.markdown(text, extensions=["tables", "sane_lists"])
        for i, block in enumerate(blocks):
            out = out.replace(f"<p>CODEBLOCK{i}XX</p>", block)
        out = out.replace("<table>", '<table class="doc-table">')
        return out

    # -- section splitting --------------------------------------------------
    @staticmethod
    def split(text, level):
        """Split text on headings of exactly `level` #'s -> [(heading, body)]."""
        marker = "#" * level + " "
        parts, head, buf = [], None, []
        in_code = False
        for line in text.splitlines():
            if line.startswith("```"):
                in_code = not in_code
            if not in_code and line.startswith(marker):
                parts.append((head, "\n".join(buf)))
                head, buf = line[len(marker):].strip(), []
            else:
                buf.append(line)
        parts.append((head, "\n".join(buf)))
        return parts

    def heading_html(self, text):
        return markdown.markdown(text).removeprefix("<p>").removesuffix("</p>")

    def render_args(self, body, level=3, prefix="arg-"):
        parts = self.split(body, level)
        out = [self.md(parts[0][1])] if parts[0][1].strip() else []
        for head, sub in parts[1:]:
            nested = self.split(sub, level + 1)
            head_html = self.heading_html(head)
            name = re.search(r"<code>(.*?)</code>", head_html)
            anchor = prefix + slug(name.group(1) if name else head)
            m = re.match(r"^(.*?</code>)\s*(?:—|--|-)\s*(.*)$", head_html)
            if m:
                summary = f'{m.group(1)} <span class="arg-desc">— {m.group(2)}</span>'
            else:
                summary = head_html
            inner = self.md(nested[0][1])
            # A first paragraph of the form *type* | *type* is the data-type line.
            inner = re.sub(r"^<p><em>(.*?)</em></p>",
                           r'<p class="arg-type">\1</p>', inner, count=1, flags=re.S)
            if len(nested) > 1:
                inner += self.render_args(sub[len(nested[0][1]):], level + 1, prefix)
            out.append(f'<details class="arg" id="{anchor}"><summary>{summary}</summary>'
                       f'<div class="arg-body">{inner}</div></details>')
        return "\n".join(out)

    def render_examples(self, body):
        parts = self.split(body, 3)
        out = [self.md(parts[0][1])] if parts[0][1].strip() else []
        for i, (head, sub) in enumerate(parts[1:]):
            anchor = "example-" + slug(head)
            open_attr = " open" if i == 0 else ""
            out.append(f'<details class="example" id="{anchor}"{open_attr}>'
                       f'<summary>{self.heading_html(head)}</summary>'
                       f'<div class="example-body">{self.md(sub)}</div></details>')
        return "\n".join(out)

    def render_history(self, body):
        parts = self.split(body, 3)
        out = [self.md(parts[0][1])] if parts[0][1].strip() else []
        for head, sub in parts[1:]:
            out.append(f'<div class="history"><h3>{self.heading_html(head)}</h3>'
                       f'{self.md(sub)}</div>')
        return "\n".join(out)

    def render_syntax(self, body, descriptions):
        code = re.search(r"```\w*\n(.*?)```", body, re.S)
        lines = [l for l in (code.group(1) if code else body).splitlines() if l.strip()]
        out = ['<div class="syntax-block">']
        for line in lines:
            key = re.sub(r"\s+", "", line.split("%")[0])
            target = descriptions.get(key)
            text = html.escape(line.strip(), quote=False)
            if target:
                out.append(f'<div class="syntax-line"><a href="#{target}"><code>{text}</code></a></div>')
            else:
                out.append(f'<div class="syntax-line"><code>{text}</code></div>')
        out.append("</div>")
        return "\n".join(out)

    @staticmethod
    def mark_descriptions(desc_html):
        """Give description paragraphs that open with a calling form an id."""
        found = {}

        def repl(m):
            code = html.unescape(m.group(1))
            key = re.sub(r"\s+", "", code)
            anchor = f"desc-{len(found) + 1}"
            found[key] = anchor
            return f'<p class="syntax-desc" id="{anchor}"><code>{m.group(1)}</code>'
        desc_html = re.sub(r"<p><code>(.*?)</code>", repl, desc_html)
        return desc_html, found

    # -- whole page -----------------------------------------------------------
    def render_body(self, page):
        sections = self.split(page.body, 2)
        rendered, toc = [], []
        intro = sections[0][1]
        if intro.strip():
            rendered.append(self.md(intro))
        descriptions = {}
        desc_html = None
        for head, body in sections[1:]:
            if head == "Description":
                desc_html, descriptions = self.mark_descriptions(self.md(body))
        for head, body in sections[1:]:
            anchor = slug(head)
            label = strip_tags(self.heading_html(head)).strip()
            toc.append((anchor, label))
            if head == "Syntax":
                inner = self.render_syntax(body, descriptions)
            elif head == "Description":
                inner = desc_html
            elif head == "Examples":
                inner = ('<p class="example-tools"><button type="button" class="toggle-all" '
                         'data-target="example">Expand all</button></p>'
                         + self.render_examples(body))
            elif head in ARG_SECTIONS:
                inner = ('<p class="example-tools"><button type="button" class="toggle-all" '
                         'data-target="arg">Expand all</button></p>'
                         + self.render_args(body, prefix={"Properties": "prop-",
                                                         "Output Arguments": "out-"}.get(head, "arg-")))
            elif head == "Version History":
                inner = self.render_history(body)
            elif head == "See Also":
                inner = f'<div class="see-also">{self.md(body)}</div>'
            else:
                inner = self.md(body)
            if head in {"See Also"} and page.stem not in ("index",):
                rendered.append(f'<section class="refsect see-also-sect" id="{anchor}">'
                                f'<h2>{self.heading_html(head)}</h2>{inner}</section>')
            else:
                rendered.append(f'<section class="refsect" id="{anchor}">'
                                f'<h2>{self.heading_html(head)}</h2>{inner}</section>')
        body_html = "\n".join(rendered)
        body_html = self.link_code(body_html, page)
        # Give remaining h3 headings ids for deep links.
        body_html = re.sub(
            r"<h3>(.*?)</h3>",
            lambda m: f'<h3 id="{slug(m.group(1))}">{m.group(1)}</h3>', body_html)
        return body_html, toc

    def option_row_ids(self, body_html):
        """Options Reference: anchor each option table row and section."""
        section = [None]

        def sect(m):
            section[0] = html.unescape(re.sub(r"<[^>]+>", "", m.group(2)))
            return m.group(0)

        def row(m):
            if section[0] and section[0].startswith("opts."):
                full = f"{section[0]}.{html.unescape(m.group(1))}"
                return f'<tr id="{slug(full)}"><td><code>{m.group(1)}</code>'
            return m.group(0)
        out = []
        for chunk in re.split(r"(<section class=\"refsect\" id=\"[^\"]*\"><h2>.*?</h2>)",
                              body_html, flags=re.S):
            m = re.match(r"<section class=\"refsect\" id=\"([^\"]*)\"><h2>(.*?)</h2>", chunk, re.S)
            if m:
                sect(m)
                out.append(chunk)
            else:
                out.append(re.sub(r"<tr>\s*<td><code>(\w+)</code>", row, chunk))
        return "".join(out)


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------
def nav_html(current):
    out = ['<nav class="toc" aria-label="Contents">',
           f'<a class="toc-home{" current" if current == "index" else ""}" href="index.html">'
           f'{TOOLBOX}</a>', '<ul>']
    for label, _, items in TOC:
        is_open = any(stem == current for stem, _ in items)
        out.append(f'<li><details{" open" if is_open else ""}><summary>{html.escape(label)}</summary><ul>')
        for stem, text in items:
            cls = ' class="current" aria-current="page"' if stem == current else ""
            out.append(f'<li><a{cls} href="{stem}.html">{html.escape(text)}</a></li>')
        out.append("</ul></details></li>")
    fl = ' class="current" aria-current="page"' if current == "FunctionList" else ""
    out.append(f'<li class="toc-flat"><a{fl} href="FunctionList.html">All Functions</a></li>')
    out.append("</ul></nav>")
    return "\n".join(out)


def breadcrumb_html(page, renderer):
    crumbs = [f'<a href="index.html">{TOOLBOX}</a>']
    section = renderer.section_of.get(page.stem)
    if page.stem == "FunctionList":
        section = None
    if section:
        landing = next(l for s, l, _ in TOC if s == section)
        if landing != page.stem:
            crumbs.append(f'<a href="{landing}.html">{html.escape(section)}</a>')
        else:
            crumbs.append(html.escape(section))
    if page.stem != "index":
        crumbs.append(f'<span>{html.escape(page.title)}</span>')
    return '<nav class="breadcrumbs" aria-label="Breadcrumb">' + \
        ' <span class="sep">›</span> '.join(crumbs) + "</nav>"


def onpage_html(toc):
    if len(toc) < 3:
        return ""
    items = "".join(f'<li><a href="#{a}">{html.escape(l)}</a></li>' for a, l in toc)
    return f'<aside class="onpage" aria-label="On this page"><h2>On this page</h2><ul>{items}</ul></aside>'


PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title} - {toolbox}</title>
<meta name="description" content="{description}">
<link rel="stylesheet" href="assets/style.css">
<script src="assets/search-index.js" defer></script>
<script src="assets/site.js" defer></script>
</head>
<body>
<header class="topbar">
  <button class="menu" type="button" aria-label="Show contents" aria-expanded="false">&#9776;</button>
  <a class="brand" href="index.html">{toolbox}</a>
  <span class="version">v{version}</span>
  <div class="search" role="search">
    <input type="search" id="search" placeholder="Search documentation" aria-label="Search documentation" autocomplete="off">
    <ul id="search-results" hidden></ul>
  </div>
</header>
<div class="layout">
{nav}
<main class="content" id="content">
{breadcrumbs}
<div class="page-head{kind}">
<h1>{title}</h1>
{purpose}
</div>
{body}
<footer class="page-foot">
<p><a href="{repo}/issues/new/choose">Report a documentation issue</a> &middot; <a href="{repo}">Source on GitHub</a></p>
<p>Copyright &copy; 2026 Mario Andr&eacute;s Mu&ntilde;oz Acosta and contributors, School of Computing and Information Systems, The University of Melbourne. Released under the PolyForm Noncommercial License 1.0.0.</p>
</footer>
</main>
{onpage}
</div>
</body>
</html>
"""


def version():
    text = (REPO_DIR / "Contents.m").read_text(encoding="utf-8")
    m = re.search(r"Version\s+([\d.]+)", text)
    return m.group(1) if m else "0"


def function_list_page(pages):
    """Generated 'All Functions' page, grouped like the TOC."""
    lines = ["# All Functions", "",
             "Every function and class in the toolbox, by category.", ""]
    for label, _, items in TOC:
        rows = [(stem, pages[stem]) for stem, _ in items
                if pages[stem].is_function or stem == "InstanceSpace"]
        if not rows:
            continue
        lines += [f"## {label}", "", "| Name | Purpose |", "|---|---|"]
        for stem, page in rows:
            lines.append(f"| [`{page.title}`]({stem}.html) | {page.purpose} |")
        lines.append("")
    path = SRC_DIR / "FunctionList.md"
    page = Page.__new__(Page)
    page.path = path
    page.stem = "FunctionList"
    page.title = "All Functions"
    page.purpose = "Every function and class in the toolbox, by category."
    page.opts_section = None
    page.body = "\n".join(lines[4:])
    page.out_name = "FunctionList.html"
    page.is_function = False
    return page


def helptoc_xml(pages):
    icon = {"Getting Started": "HelpIcon.GETTING_STARTED",
            "Release Notes": "HelpIcon.RELEASE_NOTES",
            "Topics": "HelpIcon.USER_GUIDE"}
    out = ['<?xml version="1.0" encoding="utf-8"?>',
           '<!-- Generated by doc/generate.py from its TOC. Do not edit. -->',
           '<toc version="2.0">',
           f'  <tocitem target="index.html">{html.escape(TOOLBOX)}']
    for label, landing, items in TOC:
        img = icon.get(label, "HelpIcon.FUNCTION")
        out.append(f'    <tocitem target="{items[0][0]}.html" image="{img}">{html.escape(label)}')
        for stem, text in items:
            kind = "HelpIcon.FUNCTION" if pages[stem].is_function else img
            out.append(f'      <tocitem target="{stem}.html" image="{kind}">{html.escape(text)}</tocitem>')
        out.append("    </tocitem>")
    out.append('    <tocitem target="FunctionList.html" image="HelpIcon.FUNCTION">All Functions</tocitem>')
    out += ["  </tocitem>", "</toc>", ""]
    return "\n".join(out)


def search_entry(page, body_html):
    text = re.sub(r"<pre.*?</pre>", " ", body_html, flags=re.S)
    text = re.sub(r"\s+", " ", strip_tags(text)).strip()
    heads = [strip_tags(h).strip() for h in re.findall(r"<h[23][^>]*>(.*?)</h[23]>", body_html)]
    return {"t": page.title, "u": page.out_name, "p": page.purpose,
            "h": heads, "x": text[:4000]}


def build():
    pages = load_pages()
    pages["FunctionList"] = function_list_page(pages)
    renderer = Renderer(pages)
    ver = version()
    outputs = {}
    index = []
    for stem, page in pages.items():
        body, toc = renderer.render_body(page)
        if stem == "OptionsReference":
            body = renderer.option_row_ids(body)
        purpose = f'<p class="purpose">{renderer.link_code(renderer.heading_html(page.purpose), page)}</p>' \
            if page.purpose else ""
        outputs[page.out_name] = PAGE.format(
            title=html.escape(page.title), toolbox=TOOLBOX, version=ver,
            description=html.escape(page.purpose or page.title),
            nav=nav_html(stem), breadcrumbs=breadcrumb_html(page, renderer),
            kind=" function" if page.is_function else "",
            purpose=purpose, body=body, onpage=onpage_html(toc), repo=REPO_URL)
        index.append(search_entry(page, body))
    outputs["helptoc.xml"] = helptoc_xml(pages)
    index.sort(key=lambda e: e["u"])
    outputs["assets/search-index.js"] = (
        "// Generated by doc/generate.py. Do not edit.\n"
        "window.ISA_SEARCH = " + json.dumps(index, ensure_ascii=False, separators=(",", ":")) + ";\n")
    for asset in sorted(ASSET_DIR.iterdir()):
        outputs[f"assets/{asset.name}"] = asset.read_text(encoding="utf-8")
    # GitHub Pages: serve files as-is, and a friendly 404.
    outputs[".nojekyll"] = ""
    outputs["404.html"] = outputs["index.html"].replace(
        '<div class="page-head">', '<div class="page-head"><p class="notfound">'
        'That page does not exist. Use the contents or search to find what you need.</p>', 1)
    return outputs


def main():
    check = "--check" in sys.argv[1:]
    outputs = build()
    existing = {p.relative_to(OUT_DIR).as_posix() for p in OUT_DIR.rglob("*") if p.is_file()} \
        if OUT_DIR.exists() else set()
    stale = sorted(existing - set(outputs))
    changed = sorted(name for name, text in outputs.items()
                     if not (OUT_DIR / name).is_file()
                     or (OUT_DIR / name).read_text(encoding="utf-8") != text)
    if check:
        if changed or stale:
            print("doc/html is out of date with doc/src. Run: python doc/generate.py")
            for name in changed:
                print(f"  changed: {name}")
            for name in stale:
                print(f"  stale:   {name}")
            sys.exit(1)
        print(f"doc/html is up to date ({len(outputs)} files).")
        return
    for name in stale:
        (OUT_DIR / name).unlink()
    for name, text in outputs.items():
        path = OUT_DIR / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    print(f"Wrote {len(outputs)} files to {OUT_DIR.relative_to(REPO_DIR)} "
          f"({len(changed)} changed, {len(stale)} removed).")


if __name__ == "__main__":
    main()
