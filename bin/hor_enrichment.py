#!/usr/bin/env python3
"""
hor_enrichment.py - per-chromosome read enrichment at HOR / satellite regions.

Replaces the bigWigToBedGraph -> awk -> sort-bed -> bedmap --sum chain with
direct, base-weighted sums from the bigwigs:

    region signal = sum over intervals of ( exact mean coverage x interval length )

which fixes the bedmap issue where merged equal-value bedgraph runs are counted
once regardless of length (IP and input merge differently, so their ratio was
biased). NAN-on-empty regions are handled as zero.

For each IP/control pair and each region set it reports, per chromosome:

    enrichment = IP signal in set / control signal in set

and also an "outside" set (the per-chromosome complement of the first region
set, i.e. chromosome arms) so "Within HORs" and "Outside HORs" panels come from
one run. Output: one TSV + one lollipop PNG per pair per set (same style as the
notebook figures: dashed line at 1, chromosomes sorted by score).

Usage:
  hor_enrichment.py \
      --pairs pairs.tsv \
      --regions HOR=hor.sorted.bed [--regions HSat2=hsat2.bed --regions HSat3=hsat3.bed] \
      --chrom-sizes chm13v2.0.chrom.sizes \
      --out results/_enrichment [--no-outside] [--exclude-chroms chrM]

pairs.tsv (tab-separated, no header; '#' comments allowed):
    label <TAB> ip_bigwig <TAB> control_bigwig
e.g.
    ICRF_vs_DMSO      results/POL2_ICRF/tracks/POL2_ICRF.rpkm.bw   results/POL2_DMSO/tracks/POL2_DMSO.rpkm.bw
    RAD21_WT_vs_Input results/IP_WT/tracks/IP_WT.spike.bw          results/Input_WT/tracks/Input_WT.spike.bw

The control can be a matched input, IgG, or another IP (e.g. DMSO) - the score
is simply the ratio of the two normalised signals, so use bigwigs produced with
the SAME normalisation (both RPKM, or both spike-scaled).

Backward compatibility: --compiled-pairs takes the old bedmap outputs
(compiled_*_HOR2.bed, 'chrom start end score' with optional '|' separators)
and ratios them the same way - useful to cross-check old scores against the
new base-weighted ones. Unlike the legacy row-by-row division, rows are
joined by chromosome name, so a missing chromosome cannot mispair the rest.

Dependencies: pyBigWig, matplotlib (both ship with the deeptools conda env).
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

try:
    import pyBigWig
except ImportError:
    sys.exit("ERROR: pyBigWig not found (conda install -c bioconda pybigwig, "
             "or use the deeptools environment)")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    if not sizes:
        sys.exit(f"ERROR: no chromosomes read from {path}")
    return sizes


def read_bed(path):
    """chrom -> merged, sorted [start, end) intervals."""
    raw = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                continue
            if e > s:
                raw[f[0]].append((s, e))
    merged = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        out = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[chrom] = [(s, e) for s, e in out]
    return merged


def complement(intervals, sizes):
    """per-chromosome complement of `intervals`, only for chroms present in it."""
    comp = {}
    for chrom, ivs in intervals.items():
        size = sizes.get(chrom)
        if size is None:
            continue
        out, pos = [], 0
        for s, e in ivs:
            s, e = max(0, min(s, size)), max(0, min(e, size))
            if s > pos:
                out.append((pos, s))
            pos = max(pos, e)
        if pos < size:
            out.append((pos, size))
        if out:
            comp[chrom] = out
    return comp


def region_sums(bw, intervals, sizes):
    """chrom -> base-weighted signal sum over that chromosome's intervals.

    Uses stats type='sum' (exact base-weighted sum over covered bases) when the
    installed pyBigWig supports it, else falls back to exact mean x length
    (identical for gap-free bamCoverage bigwigs, which write zeros explicitly).
    """
    sums = {}
    bw_chroms = bw.chroms()
    for chrom, ivs in intervals.items():
        if chrom not in bw_chroms:
            continue
        limit = min(sizes.get(chrom, bw_chroms[chrom]), bw_chroms[chrom])
        total = 0.0
        for s, e in ivs:
            s, e = max(0, min(s, limit)), max(0, min(e, limit))
            if e <= s:
                continue
            try:
                v = bw.stats(chrom, s, e, type="sum", exact=True)[0]
            except RuntimeError:
                m = bw.stats(chrom, s, e, type="mean", exact=True)[0]
                v = m * (e - s) if m is not None else None
            if v is not None:
                total += v
        sums[chrom] = total
    return sums


def chrom_sort_key(name):
    base = name[3:] if name.lower().startswith("chr") else name
    try:
        return (0, int(base), "")
    except ValueError:
        return (1, 0, base)


def lollipop(scores, title, out_png, color):
    items = [(c, v) for c, v in scores.items() if v is not None]
    if not items:
        return False
    items.sort(key=lambda kv: kv[1], reverse=True)
    chroms = [c for c, _ in items][::-1]          # biggest at top
    vals = [v for _, v in items][::-1]
    fig_h = max(2.5, 0.28 * len(items) + 1.2)
    fig, ax = plt.subplots(figsize=(4.2, fig_h))
    y = range(len(chroms))
    ax.hlines(y, 0, vals, color=color, alpha=0.45, linewidth=5)
    ax.plot(vals, list(y), "o", color=color, markersize=4.5)
    ax.axvline(1, linestyle="--", color="grey", linewidth=0.8)
    ax.set_yticks(list(y))
    ax.set_yticklabels(chroms, fontsize=6)
    ax.set_xlabel("enrichment score", fontsize=7)
    ax.set_ylabel("chromosome", fontsize=7)
    ax.set_title(title, fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
    xmax = max(2.0, max(vals) * 1.15)
    ax.set_xlim(0, xmax)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs",
                    help="TSV: label, IP bigwig, control bigwig")
    ap.add_argument("--compiled-pairs",
                    help="TSV: label, tx compiled bed, ctrl compiled bed "
                         "(old bedmap 'chrom start end score' outputs)")
    ap.add_argument("--regions", action="append",
                    metavar="NAME=BED", help="region set; repeatable")
    ap.add_argument("--chrom-sizes",
                    help="required with --pairs")
    ap.add_argument("--out", default="results/_enrichment")
    ap.add_argument("--no-outside", action="store_true",
                    help="skip the outside/chrArm complement of the first region set")
    ap.add_argument("--outside-from", default=None,
                    help="region set to complement for 'outside' (default: first)")
    ap.add_argument("--exclude-chroms", default="chrM",
                    help="comma list dropped from every analysis (default chrM)")
    ap.add_argument("--colors", default="#c05555,#7bae7b,#6b8fc9,#e0a54c",
                    help="comma list of plot colours, cycled per region set")
    args = ap.parse_args()

    if not args.pairs and not args.compiled_pairs:
        ap.error("need --pairs (bigwigs) and/or --compiled-pairs (old bedmap files)")
    if args.pairs and not (args.regions and args.chrom_sizes):
        ap.error("--pairs needs --regions and --chrom-sizes")

    excluded = {c for c in args.exclude_chroms.split(",") if c}
    os.makedirs(args.out, exist_ok=True)
    colors = args.colors.split(",")
    tsv_path = os.path.join(args.out, "enrichment.tsv")
    out_fh = open(tsv_path, "w", newline="")
    w = csv.writer(out_fh, delimiter="\t")
    w.writerow(["pair", "region_set", "chrom", "ip_signal",
                "control_signal", "enrichment"])

    # ---------------------------------------------- old compiled-bed mode
    if args.compiled_pairs:
        def read_compiled(path):
            vals = {}
            with open(path) as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    f = [x for x in line.replace("|", "\t").split() if x]
                    if len(f) < 4:
                        continue
                    try:
                        vals[f[0]] = float(f[3])
                    except ValueError:
                        vals[f[0]] = 0.0        # bedmap NAN and friends
            return vals

        with open(args.compiled_pairs) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    sys.exit(f"ERROR: bad compiled-pairs line: {line!r}")
                label, tx_path, ctl_path = f[0], f[1], f[2]
                tx, ctl = read_compiled(tx_path), read_compiled(ctl_path)
                scores = {}
                for chrom in sorted(set(tx) | set(ctl), key=chrom_sort_key):
                    if chrom in excluded:
                        continue
                    a, b = tx.get(chrom, 0.0), ctl.get(chrom, 0.0)
                    score = (a / b) if b > 0 else None
                    scores[chrom] = score
                    w.writerow([label, "compiled", chrom, f"{a:.4f}", f"{b:.4f}",
                                f"{score:.4f}" if score is not None else "NA"])
                png = os.path.join(args.out, f"{label}.compiled.png")
                lollipop(scores, f"Read enrichment {label} (compiled)",
                         png, colors[0])
                print(f"{label} / compiled: plot {png}")
        if not args.pairs:
            out_fh.close()
            print(f"table: {tsv_path}")
            return

    sizes = read_chrom_sizes(args.chrom_sizes)

    sets = []
    for spec in args.regions:
        if "=" not in spec:
            sys.exit(f"ERROR: --regions needs NAME=BED, got '{spec}'")
        name, bed = spec.split("=", 1)
        ivs = read_bed(bed)
        ivs = {c: v for c, v in ivs.items() if c not in excluded}
        if not ivs:
            sys.exit(f"ERROR: no usable intervals in {bed}")
        sets.append((name, ivs))

    if not args.no_outside:
        base_name = args.outside_from or sets[0][0]
        base = dict(sets)[base_name] if base_name in dict(sets) else None
        if base is None:
            sys.exit(f"ERROR: --outside-from '{base_name}' is not a loaded region set")
        sets.append((f"outside_{base_name}", complement(base, sizes)))

    pairs = []
    with open(args.pairs) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                sys.exit(f"ERROR: bad pairs line (need 3 tab-separated fields): {line!r}")
            for p in (f[1], f[2]):
                if not os.path.isfile(p):
                    sys.exit(f"ERROR: bigwig not found: {p}")
            pairs.append((f[0], f[1], f[2]))
    if not pairs:
        sys.exit("ERROR: no pairs read")

    for label, ip_path, ctl_path in pairs:
        bw_ip = pyBigWig.open(ip_path)
        bw_ctl = pyBigWig.open(ctl_path)
        for si, (set_name, ivs) in enumerate(sets):
            ip_sums = region_sums(bw_ip, ivs, sizes)
            ctl_sums = region_sums(bw_ctl, ivs, sizes)
            scores = {}
            g_ip = g_ctl = 0.0
            for chrom in sorted(ivs, key=chrom_sort_key):
                a = ip_sums.get(chrom)
                b = ctl_sums.get(chrom)
                if a is None and b is None:
                    continue                      # chrom absent from both bigwigs
                a, b = a or 0.0, b or 0.0
                g_ip += a
                g_ctl += b
                score = (a / b) if b > 0 else None
                scores[chrom] = score
                w.writerow([label, set_name, chrom, f"{a:.4f}", f"{b:.4f}",
                            f"{score:.4f}" if score is not None else "NA"])
            w.writerow([label, set_name, "ALL", f"{g_ip:.4f}", f"{g_ctl:.4f}",
                        f"{g_ip/g_ctl:.4f}" if g_ctl > 0 else "NA"])
            png = os.path.join(args.out, f"{label}.{set_name}.png")
            ok = lollipop(scores, f"Read enrichment {label} ({set_name})",
                          png, colors[si % len(colors)])
            print(f"{label} / {set_name}: "
                  f"{'plot ' + png if ok else 'no scores to plot'}")
        bw_ip.close()
        bw_ctl.close()

    out_fh.close()
    print(f"table: {tsv_path}")


if __name__ == "__main__":
    main()
