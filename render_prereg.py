#!/usr/bin/env python3
"""Render a preregistration markdown file to a clean, readable PDF.

Handles: # / ## / ### headings, **bold** (rendered plain), `code`
(rendered plain), - bullet lists, plain paragraphs, --- rules.
Long paragraphs flow across pages automatically.
Built for Abby's preregistration documents (Oct 2026).
Usage: python3 render_prereg.py input.md output.pdf
"""
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

NAVY, GOLD = "#1b2a4a", "#b98a1e"
PAGE_W, PAGE_H = 8.5, 11
LEFT = 0.85 / PAGE_W          # normalized x of left margin
TOP = (PAGE_H - 0.95) / PAGE_H
BOTTOM = 0.78 / PAGE_H
LINE_H = 0.0245               # normalized line height at 10pt


def parse(md_text):
    blocks, para = [], []

    def flush():
        if para:
            blocks.append(("para", " ".join(para)))
            para.clear()

    for raw in md_text.splitlines():
        line = raw.rstrip()
        if not line.strip():
            flush()
            continue
        if line.strip() == "---":
            flush()
            blocks.append(("rule", ""))
            continue
        m = re.match(r"^(#{1,3})\s+(.*)$", line)
        if m:
            flush()
            blocks.append((f"h{len(m.group(1))}", m.group(2).strip()))
            continue
        if re.match(r"^\s*[-*]\s+", line):
            flush()
            blocks.append(("bullet", re.sub(r"^\s*[-*]\s+", "", line).strip()))
            continue
        para.append(line.strip())
    flush()
    return blocks


def plain(text):
    text = re.sub(r"\*\*(.+?)\*\*", r"\1", text)
    text = re.sub(r"`(.+?)`", r"\1", text)
    return text


def wrap(text, width):
    words, cur, lines = text.split(" "), [], []
    for w in words:
        if len(" ".join(cur + [w])) > width:
            lines.append(" ".join(cur))
            cur = [w]
        else:
            cur.append(w)
    if cur:
        lines.append(" ".join(cur))
    return lines


class Doc:
    def __init__(self, pdf, footer):
        self.pdf, self.footer = pdf, footer
        self.fig = None
        self.new_page()

    def new_page(self):
        if self.fig is not None:
            self.fig.text(0.5, 0.38 / PAGE_H, self.footer, fontsize=8,
                         ha="center", color="#888888")
            self.pdf.savefig(self.fig)
            plt.close(self.fig)
        self.fig = plt.figure(figsize=(PAGE_W, PAGE_H))
        self.fig.patch.set_facecolor("white")
        self.y = TOP

    def ensure(self, lines_needed=2, fs=10.0):
        if self.y - lines_needed * LINE_H * (fs / 10.0) < BOTTOM:
            self.new_page()

    def line(self, x, s, fs=10.0, weight="normal", color="#222222"):
        if self.y < BOTTOM:
            self.new_page()
        self.fig.text(x, self.y, s, fontsize=fs, va="top", ha="left",
                      weight=weight, color=color)
        self.y -= LINE_H * (fs / 10.0)

    def para(self, text, fs=10.0, x=LEFT, width=92, gap=0.10):
        for ln in wrap(plain(text), int(width * (10.0 / fs))):
            self.line(x, ln, fs)
        self.y -= gap / PAGE_H

    def finish(self):
        self.fig.text(0.5, 0.38 / PAGE_H, self.footer, fontsize=8,
                     ha="center", color="#888888")
        self.pdf.savefig(self.fig)
        plt.close(self.fig)


def render(md_path, pdf_path, footer="Abby Davis · preregistration · 2026-10-03"):
    with open(md_path) as f:
        blocks = parse(f.read())
    with PdfPages(pdf_path) as pdf:
        d = Doc(pdf, footer)
        for btype, text in blocks:
            if btype == "rule":
                d.ensure(2)
                ax = d.fig.add_axes([LEFT, d.y - 0.012, 1 - 2 * LEFT, 0.004])
                ax.axis("off")
                ax.add_patch(plt.Rectangle((0, 0), 1, 1, color=GOLD))
                d.y -= 0.26 / PAGE_H
            elif btype == "h1":
                d.ensure(4, 17)
                for ln in wrap(plain(text), 52):
                    d.line(LEFT, ln, 17, weight="bold", color=NAVY)
                ax = d.fig.add_axes([LEFT, d.y - 0.010, 1 - 2 * LEFT, 0.004])
                ax.axis("off")
                ax.add_patch(plt.Rectangle((0, 0), 1, 1, color=GOLD))
                d.y -= 0.20 / PAGE_H
            elif btype in ("h2", "h3"):
                fs = 13 if btype == "h2" else 11.5
                d.ensure(3, fs)
                for ln in wrap(plain(text), int(92 * (10.0 / fs))):
                    d.line(LEFT, ln, fs, weight="bold", color=NAVY)
                d.y -= 0.06 / PAGE_H
            elif btype == "bullet":
                d.ensure(2)
                d.line(LEFT - 0.022, "•", 10.0, weight="bold", color=GOLD)
                # pull the bullet line back up: draw text at same y
                d.y += LINE_H
                d.para(text, x=LEFT + 0.004, width=90, gap=0.08)
            else:
                d.ensure(2)
                d.para(text)
        d.finish()
    print("wrote", pdf_path)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: render_prereg.py input.md output.pdf")
        sys.exit(2)
    render(sys.argv[1], sys.argv[2])
