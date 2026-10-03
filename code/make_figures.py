#!/usr/bin/env python3
"""make_figures.py — prereg v2.0 §11 figures from data/figure JSONs.

(1) zeta event map (t vs z, type-encoded); (2) fold/landing fits with
residuals; (3) N_off(t) for all cases; (4) artifact gallery.
"""
import json
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
V2DIR = os.path.dirname(HERE)
FIGDIR = os.path.join(V2DIR, 'figures')
DATADIR = os.path.join(V2DIR, 'data')


def load(name):
    with open(os.path.join(FIGDIR, name)) as f:
        return json.load(f)


def fig1():
    evs = load('fig1_event_map.json')
    fig, ax = plt.subplots(figsize=(10, 6))
    marks = {'R': ('o', 'tab:blue'), 'C': ('^', 'tab:red'),
             'X': ('x', 'tab:gray')}
    for e in evs:
        m, c = marks.get(e['type'], ('s', 'black'))
        ax.errorbar(e['t_star'], e['z_re'], xerr=e['sigma'], fmt=m,
                    color=c, label=e['type'], capsize=2, ms=6)
    # de-duplicate legend
    handles, labels = ax.get_legend_handles_labels()
    seen = {}
    for h, l in zip(handles, labels):
        seen[l] = h
    ax.legend(seen.values(), seen.keys(), title='event type')
    ax.set_xlabel('t* (fitted event time)')
    ax.set_ylabel('Re(z*)')
    ax.set_title('Fig 1 — zeta validated event map (prereg v2.0 §11)\n'
                 'empty map = reported plainly, never rescued')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig1.png'), dpi=120)
    plt.close(fig)
    print('fig1.png (%d events)' % len(evs))


def fig2():
    fits = load('fig2_fits.json')
    fits = [f for f in fits if f.get('fit_pts')]
    n = min(len(fits), 12)
    if n == 0:
        print('fig2: no fits with points; skipped')
        return
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3.2 * rows),
                             squeeze=False)
    for ax, f in zip(axes.flat, fits[:n]):
        pts = f['fit_pts']
        ts = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.scatter(ts, ys, s=10, color='black', zorder=3)
        # fit line: y = C^2 (t - t*)
        import numpy as np
        tg = np.linspace(min(ts), max(ts), 100)
        C2 = (f['C'] or 0) ** 2
        ax.plot(tg, C2 * (tg - f['t_star']), color='tab:red', lw=1.5)
        ax.axvline(f['t_star'], color='tab:red', ls='--', lw=1)
        ax.set_title('%s %s t*=%.4f R²=%.5f'
                     % (f['case'], f['type'], f['t_star'], f['R2'] or 0),
                     fontsize=9)
        ax.set_xlabel('t')
        ax.set_ylabel('g² / im²')
        ax.grid(alpha=0.3)
    for ax in axes.flat[n:]:
        ax.axis('off')
    fig.suptitle('Fig 2 — fold/landing fits (prereg v2.0 §11)')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig2.png'), dpi=120)
    plt.close(fig)
    print('fig2.png (%d fits shown of %d)' % (n, len(fits)))


def fig3():
    noff = load('fig3_noff.json')
    fig, ax = plt.subplots(figsize=(10, 6))
    for name, curve in sorted(noff.items()):
        ts = [p[0] for p in curve]
        ns = [p[1] for p in curve]
        ax.step(ts, ns, where='post', label=name)
    ax.set_xlabel('t')
    ax.set_ylabel('N_off(t)  (|Im z| > 1e-6)')
    ax.set_title('Fig 3 — off-axis zero counts (prereg v2.0 §11)')
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig3.png'), dpi=120)
    plt.close(fig)
    print('fig3.png')


def fig4():
    gal = load('fig4_artifacts.json')
    fig, axes = plt.subplots(1, 2, figsize=(13, 5),
                             gridspec_kw={'width_ratios': [1, 2]})
    # left: verdict counts per case from events.csv
    import csv
    counts = {}
    with open(os.path.join(DATADIR, 'events.csv')) as f:
        for row in csv.DictReader(f):
            key = (row['case'], row['verdict'])
            counts[key] = counts.get(key, 0) + 1
    cases = sorted({c for c, _ in counts})
    verdicts = ['genuine', 'suspect', 'artifact', 'logged']
    x = range(len(cases))
    bottom = [0] * len(cases)
    for v in verdicts:
        vals = [counts.get((c, v), 0) for c in cases]
        axes[0].bar(x, vals, bottom=bottom, label=v)
        bottom = [b + a for b, a in zip(bottom, vals)]
    axes[0].set_xticks(list(x))
    axes[0].set_xticklabels(cases, rotation=20)
    axes[0].set_ylabel('detector triggers')
    axes[0].set_title('Fig 4a — trigger verdicts per case')
    axes[0].legend()
    # right: artifact/suspect reasons table
    axes[1].axis('off')
    rows = [[g.get('case', ''), g.get('type', ''), g.get('t_detect', ''),
             str(g.get('artifact_reason', ''))[:90]] for g in gal[:18]]
    if rows:
        tbl = axes[1].table(cellText=rows,
                            colLabels=['case', 'type', 't_detect', 'reason'],
                            loc='center')
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(7)
    axes[1].set_title('Fig 4b — artifact gallery (triggers failing G-a/b/c)')
    fig.suptitle('Fig 4 — artifact gallery (prereg v2.0 §11)')
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, 'fig4.png'), dpi=120)
    plt.close(fig)
    print('fig4.png (%d artifact/suspect rows)' % len(gal))


def main():
    fig1()
    fig2()
    fig3()
    fig4()


if __name__ == '__main__':
    main()
