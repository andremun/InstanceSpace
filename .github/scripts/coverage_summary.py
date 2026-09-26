#!/usr/bin/env python3
"""Print a Markdown per-file line-coverage table from a Cobertura XML report.

Usage: python3 coverage_summary.py coverage.xml

Each row lists the file, covered/total executable lines, the percentage,
and the uncovered line numbers collapsed into ranges. Exits 0 even when the
report is missing, so a failed test run still reaches later CI steps.
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def ranges(nums):
    """Collapse a sorted list of ints into 'a-b' range strings."""
    out, start, prev = [], None, None
    for n in nums:
        if start is None:
            start = prev = n
        elif n == prev + 1:
            prev = n
        else:
            out.append(f"{start}" if start == prev else f"{start}-{prev}")
            start = prev = n
    if start is not None:
        out.append(f"{start}" if start == prev else f"{start}-{prev}")
    return out


def main(path):
    report = Path(path)
    if not report.is_file():
        print(f"No coverage report at {report}.")
        return
    files = {}
    for cls in ET.parse(report).getroot().iter("class"):
        name = cls.get("filename")
        lines = files.setdefault(name, {})
        for line in cls.iter("line"):
            num = int(line.get("number"))
            lines[num] = lines.get(num, 0) + int(line.get("hits"))
    total_hit = total = 0
    rows = []
    for name in sorted(files):
        lines = files[name]
        hit = sum(1 for h in lines.values() if h > 0)
        missed = sorted(n for n, h in lines.items() if h == 0)
        total_hit += hit
        total += len(lines)
        pct = 100.0 * hit / len(lines) if lines else 100.0
        rows.append((name, hit, len(lines), pct, ", ".join(ranges(missed))))
    overall = 100.0 * total_hit / total if total else 0.0
    print(f"## Line coverage: {overall:.1f}% ({total_hit}/{total})\n")
    print("| File | Covered | % | Uncovered lines |")
    print("|---|---|---|---|")
    for name, hit, n, pct, missed in rows:
        print(f"| `{name}` | {hit}/{n} | {pct:.1f} | {missed} |")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "coverage.xml")
