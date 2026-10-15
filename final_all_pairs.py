# PRA2003 - all six normal/variant pairs.
# Output: counts per event, overall asymmetry, asymmetry vs P (with quantification), momentum spectrum dN/dP.
# Method (sub-sampling): the data are split into 10 similar-sized sub-samples (the 10 files). For every quantity the
# (weighted) average over the sub-samples is the result and the standard deviation (ddof=1) of the sub-sample
# results is the statistical uncertainty.
# Put output-Set1.txt ... output-Set10.txt in the same directory as this script, then run: python final_all_pairs.py
# Also writes: asymmetry_vs_P.png, spectrum_dNdP.png and results_Q3.md (tables to paste into the README).

import numpy as np
import matplotlib.pyplot as plt
from bisect import bisect_right
from multiprocessing import Pool

try:
    from scipy.stats import chi2 as chi2_dist
except ImportError:                      # p-values are skipped if scipy is not installed
    chi2_dist = None

FILE_PATTERN = "output-Set{}.txt"
FILE_NUMBERS = range(1, 11)              # Set1..Set10 (range(1, 11) = 10 files). Use range(0, 11) if you also have Set0
PAIRS = [(211, -211), (321, -321), (2212, -2212), (3122, -3122), (3312, -3312), (3334, -3334)]  # (normal, variant)
NAMES = {211: "Carbon monoxide", -211: "Carbon-13 monoxide",
         321: "Nitric oxide", -321: "Ionised NO",
         2212: "Water", -2212: "Heavy water (D2O)",
         3122: "Methane (CH4)", -3122: "Methyl ion (CH3-)",
         3312: "Ethylene (C2H4)", -3312: "Ionised ethylene (C2H3-)",
         3334: "Ozone (O3)", -3334: "Superoxide anion (O2-)"}

# True: counts per event are divided by ALL events (including empty ones). False: only non-empty events.
INCLUDE_EMPTY_EVENTS = False

# Fine bin edges (variable width). The last bin is an open tail: everything above the last finite edge.
BINS_FINE = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, np.inf])
N_FINE = len(BINS_FINE) - 1
# Rare species (ethylene, ozone) have too few particles per sub-sample for the fine bins, so adjacent fine bins are
# merged into these coarser bins: 0-1, 1-2, 2-3, 3-6, >6. Indices refer to BINS_FINE.
COARSE_START = [0, 2, 4, 6, 8]
COARSE_PAIRS = {3312, 3334}
P_LABEL = "P ($10^{-23}$ kg m/s)"
BINS_LIST = BINS_FINE.tolist()


def geometry(bins):
    """Plotting geometry: the open tail bin is drawn as wide as the last finite bin."""
    edges = bins.copy()
    edges[-1] = bins[-2] + (bins[-2] - bins[-3])
    return edges, (edges[:-1] + edges[1:]) / 2, np.diff(edges)


def analyse_file(filename):
    """Returns (events used, total counts per ID, counts per ID per fine momentum bin)."""
    n_events = 0
    totals = {}
    binned = {}
    with open(filename, "r") as infile:
        while True:
            header = infile.readline().split()
            if not header:
                break
            n_particles = int(header[1])
            if n_particles == 0 and not INCLUDE_EMPTY_EVENTS:
                continue
            n_events += 1
            for _ in range(n_particles):
                line = infile.readline().split()
                code = int(line[3])
                totals[code] = totals.get(code, 0) + 1
                x, y, z = float(line[0]), float(line[1]), float(line[2])
                b = bisect_right(BINS_LIST, (x * x + y * y + z * z) ** 0.5) - 1
                if 0 <= b < N_FINE:
                    if code not in binned:
                        binned[code] = np.zeros(N_FINE)
                    binned[code][b] += 1
    return n_events, totals, binned


def asymmetry(n_normal, n_variant):
    """(N_normal - N_variant)/(N_normal + N_variant); nan where the total is 0. The bin width cancels."""
    a, b = np.asarray(n_normal, float), np.asarray(n_variant, float)
    total = a + b
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total > 0, (a - b) / total, np.nan)


def weighted_mean_std(x, w, axis=0):
    """Weighted mean over sub-samples and sample std (ddof=1) of the sub-samples; nan / zero-weight entries ignored."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    valid = np.isfinite(x) & (w > 0)
    xs = np.where(valid, x, 0.0)
    ws = np.where(valid, w, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.sum(ws * xs, axis=axis) / np.sum(ws, axis=axis)
        std = np.nanstd(np.where(valid, x, np.nan), axis=axis, ddof=1)
    return mean, std


def p_value(chi2, ndof):
    return chi2_dist.sf(chi2, ndof) if (chi2_dist is not None and ndof > 0) else np.nan


def merge(counts, normal):
    """counts: [file][fine bin]. Returns counts in the bins used for this pair, plus the bin edges."""
    if normal in COARSE_PAIRS:
        bins = np.append(BINS_FINE[COARSE_START], np.inf)
        return np.add.reduceat(counts, COARSE_START, axis=1), bins
    return counts, BINS_FINE


def main():
    files = [FILE_PATTERN.format(n) for n in FILE_NUMBERS]
    n_files = len(files)
    print("Reading", n_files, "files in parallel ...")
    with Pool(min(n_files, 10)) as pool:                       # one process per file
        read = pool.map(analyse_file, files)

    n_events_all = np.array([r[0] for r in read], float)
    n_pairs = len(PAIRS)
    total_counts = np.zeros((n_pairs, n_files))
    counts_per_event = np.zeros((n_pairs, 2, n_files))
    total_asym = np.zeros((n_pairs, n_files))
    fine_normal = np.zeros((n_pairs, n_files, N_FINE))
    fine_variant = np.zeros((n_pairs, n_files, N_FINE))
    zeros = np.zeros(N_FINE)
    for f, (n_events, totals, binned) in enumerate(read):
        for i, (normal, variant) in enumerate(PAIRS):
            total_counts[i, f] = totals.get(normal, 0) + totals.get(variant, 0)
            counts_per_event[i, 0, f] = totals.get(normal, 0) / n_events
            counts_per_event[i, 1, f] = totals.get(variant, 0) / n_events
            total_asym[i, f] = asymmetry(totals.get(normal, 0), totals.get(variant, 0))
            fine_normal[i, f] = binned.get(normal, zeros)
            fine_variant[i, f] = binned.get(variant, zeros)

    print("\nEvents per sub-sample:", n_events_all.astype(int))
    print("(max/min = %.3f; close to 1 means the sub-samples are similar in size)"
          % (n_events_all.max() / n_events_all.min()))
    print("Empty events %s in the denominator.\n" % ("INCLUDED" if INCLUDE_EMPTY_EVENTS else "EXCLUDED"))

    md = ["# Question 3 results (generated by final_all_pairs.py)\n"]

    # ---------------- overall: counts per event and asymmetry ----------------
    overall = []
    md.append("## Overall asymmetry per pair\n")
    md.append("| Pair | Asymmetry (+- std of sub-samples) | Significance |")
    md.append("|---|---|---|")
    for i, (normal, variant) in enumerate(PAIRS):
        for j, code in enumerate((normal, variant)):
            m, s = weighted_mean_std(counts_per_event[i, j], n_events_all)
            print(f"ID {code} ({NAMES[code]}): {m:.5f} +- {s:.5f} counts per event")
        m, s = weighted_mean_std(total_asym[i], total_counts[i])
        overall.append((m, s))
        print(f"Overall asymmetry {normal}/{variant}: {m:+.5f} +- {s:.5f}  ({abs(m) / s:.1f} sigma from zero)\n")
        md.append(f"| {NAMES[normal]} vs {NAMES[variant]} | {m:+.5f} +- {s:.5f} | {abs(m) / s:.1f} sigma |")

    # ---------------- asymmetry vs P ----------------
    res = []
    md.append("\n## Asymmetry as a function of momentum\n")
    for i, (normal, variant) in enumerate(PAIRS):
        bn, bins = merge(fine_normal[i], normal)
        bv, _ = merge(fine_variant[i], normal)
        bc = bn + bv
        asym = asymmetry(bn, bv)
        mean, std = weighted_mean_std(asym, bc, axis=0)
        edges, centers, widths = geometry(bins)
        nb = len(bins) - 1
        res.append(dict(bins=bins, edges=edges, centers=centers, widths=widths, mean=mean, std=std,
                        bn=bn, bv=bv, bc=bc))

        print(f"\nPair {normal}/{variant}: smallest N + N' per bin (min over sub-samples):",
              bc.min(axis=0).astype(int))
        print("  P range          asymmetry +- std        |A|/std")
        md.append(f"### {NAMES[normal]} vs {NAMES[variant]} ({normal}/{variant})\n")
        md.append("| P range (1e-23 kg m/s) | Asymmetry +- std | \\|A\\|/std |")
        md.append("|---|---|---|")
        for b in range(nb):
            sig = abs(mean[b]) / std[b] if std[b] > 0 else np.nan
            print(f"  {bins[b]:>4g} - {bins[b + 1]:<5g}   {mean[b]:+.5f} +- {std[b]:.5f}   {sig:5.1f}")
            md.append(f"| {bins[b]:g} - {bins[b + 1]:g} | {mean[b]:+.5f} +- {std[b]:.5f} | {sig:.1f} |")

        ok = np.isfinite(mean) & (std > 0)
        n = int(ok.sum())
        a, s = mean[ok], std[ok]
        a_all = overall[i][0]
        chi2_zero = np.sum((a / s) ** 2)
        chi2_const = np.sum(((a - a_all) / s) ** 2)
        fin = ok.copy()
        fin[-1] = False                                           # drop the open tail bin from the line fit
        coef, cov = np.polyfit(centers[fin], mean[fin], 1, w=1 / std[fin], cov="unscaled")
        slope, slope_err = coef[0], np.sqrt(cov[0, 0])
        res[-1].update(chi2_zero=chi2_zero, n=n)
        sig_bins = np.where(ok & (np.abs(mean / np.where(std > 0, std, np.nan)) > 3))[0]
        sig_txt = ", ".join(f"{bins[b]:g}-{bins[b + 1]:g}" for b in sig_bins) if len(sig_bins) else "none"
        print(f"  chi2 vs zero asymmetry : {chi2_zero:7.2f} / {n}   p = {p_value(chi2_zero, n):.3g}")
        print(f"  chi2 vs constant ({a_all:+.4f}) : {chi2_const:7.2f} / {n}   p = {p_value(chi2_const, n):.3g}")
        print(f"  slope dA/dP : {slope:+.5f} +- {slope_err:.5f} ({abs(slope) / slope_err:.1f} sigma from 0)")
        print("  bins with |A|/std > 3 :", sig_txt)
        md.append("")
        md.append(f"- chi2 vs zero asymmetry: {chi2_zero:.2f} / {n}, p = {p_value(chi2_zero, n):.3g}")
        md.append(f"- chi2 vs constant overall asymmetry ({a_all:+.4f}): {chi2_const:.2f} / {n}, "
                  f"p = {p_value(chi2_const, n):.3g}")
        md.append(f"- linear slope dA/dP = {slope:+.5f} +- {slope_err:.5f} per 1e-23 kg m/s "
                  f"({abs(slope) / slope_err:.1f} sigma from 0)")
        md.append(f"- bins with |A|/std > 3: {sig_txt}\n")

    # ---------------- homogeneity check ----------------
    print("\n" + "=" * 70 + "\nSub-sample homogeneity check (counts N + N' per event, per bin)\n" + "=" * 70)
    for i, (normal, variant) in enumerate(PAIRS):
        bc = res[i]["bc"]
        rel = bc / n_events_all[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            z = (rel - rel.mean(axis=0)) / rel.std(axis=0, ddof=1)
        flagged = [number for f, number in enumerate(FILE_NUMBERS) if np.nanmax(np.abs(z[f])) > 2.5]
        print(f"Pair {normal}/{variant}: files with |z| > 2.5 in some bin:", flagged if flagged else "none")

    # ---------------- plot: asymmetry vs P (2 x 3 grid) ----------------
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for i, (normal, variant) in enumerate(PAIRS):
        ax = axes.flat[i]
        r = res[i]
        ax.errorbar(r["centers"], r["mean"], yerr=r["std"], xerr=r["widths"] / 2, fmt="o", color="tab:blue",
                    capsize=3, label="Weighted mean of sub-sample asymmetries\n(vertical bar = std of sub-samples,\n"
                                     "horizontal bar = bin width)")
        ax.axhline(0, color="gray", ls="--", lw=0.8)
        m, sd = overall[i]
        ax.axhspan(m - sd, m + sd, color="tab:orange", alpha=0.25, label="Overall asymmetry $\\pm$ std")
        ax.text(0.02, 0.97, f"$\\chi^2$ vs 0: {r['chi2_zero']:.1f} / {r['n']}", transform=ax.transAxes,
                ha="left", va="top", fontsize=9)
        ax.set_xticks(r["edges"][:-1])
        ax.set_xticklabels([f"{e:g}" for e in r["bins"][:-1]])
        ax.set_xlim(r["edges"][0] - 0.25, r["edges"][-1] + 0.25)
        top = np.nanmax(np.abs(r["mean"]) + np.nan_to_num(r["std"]))
        ax.set_ylim(-1.2 * top, 1.2 * top)
        ax.text(0.98, 0.02, f"last bin: P > {r['bins'][-2]:g}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color="gray")
        ax.set_xlabel(P_LABEL)
        ax.set_ylabel(f"(N({normal}) - N({variant})) / (N({normal}) + N({variant}))", fontsize=8)
        ax.set_title(f"{NAMES[normal]} / {NAMES[variant]}", fontsize=10)
        ax.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    plt.tight_layout()
    plt.savefig("asymmetry_vs_P.png", dpi=150, bbox_inches="tight")

    # ---------------- plot: momentum spectrum dN/dP ----------------
    fig2, axes2 = plt.subplots(2, 3, figsize=(16, 8))
    for i, (normal, variant) in enumerate(PAIRS):
        ax = axes2.flat[i]
        r = res[i]
        for counts, code, colour in ((r["bn"], normal, "tab:blue"), (r["bv"], variant, "tab:red")):
            dens = counts / n_events_all[:, None] / r["widths"]
            w = np.broadcast_to(n_events_all[:, None], dens.shape)
            m, s = weighted_mean_std(dens, w, axis=0)
            ax.errorbar(r["centers"][:-1], m[:-1], yerr=s[:-1], xerr=r["widths"][:-1] / 2, fmt="o",
                        color=colour, capsize=3, label=f"ID {code} ({NAMES[code]})")
        ax.set_yscale("log")
        ax.set_xlabel(P_LABEL)
        ax.set_ylabel("dN/dP (counts per event per unit P)")
        ax.set_title(f"{NAMES[normal]} / {NAMES[variant]}", fontsize=10)
        ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig("spectrum_dNdP.png", dpi=150)

    with open("results_Q3.md", "w") as out:
        out.write("\n".join(md) + "\n")
    print("\nSaved asymmetry_vs_P.png, spectrum_dNdP.png and results_Q3.md")
    plt.show()


if __name__ == "__main__":                 # required for multiprocessing on Windows / macOS
    main()
