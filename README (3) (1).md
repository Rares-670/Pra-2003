# PRA2003 - Simulating molecular emissions in a combustion reaction

Rares Weintraub - i6438440

Answer the following questions:

1. What are the average counts of each molecular species and their statistical uncertainties?
2. Is there any asymmetry between the normal and the variant molecule?
3. Is there any asymmetry as a function of their momentum?

## Installation and usage

**Requirements:** Python 3 with `numpy` and `matplotlib`. `scipy` is optional (only used to print p-values for the χ² tests; the script runs without it).

```
pip install numpy matplotlib scipy
```

1. Download the data files from the surfdrive link (https://surfdrive.surf.nl/index.php/s/7udCnWTk4yMUASD) and place them in the same folder as `final_all_pairs.py`, named `output-Set1.txt` to `output-Set10.txt`.
2. Open a terminal in that folder and run: `python final_all_pairs.py`
3. The script reads the 10 files in parallel (one process per file) and:
   - prints the average counts per event of each species with their uncertainties and the overall asymmetry of each pair (Questions 1 and 2),
   - prints the asymmetry in every momentum bin, a χ² test against zero and against a constant asymmetry, a linear-slope fit, and a sub-sample homogeneity check (Question 3),
   - saves the plots `asymmetry_vs_P.png` and `spectrum_dNdP.png`, and the per-bin tables in `results_Q3.md`.

Runtime: about 3 minutes on my machine (each file is ~770 MB).

## Context and scientific question

The simulation of a combustion reaction is analysed, where a certain number of molecules is produced per event (one reaction). Each data file contains 500,000 events. For every event the file gives the event number and the number of molecules produced, followed by one line per molecule with its three momentum components (Px, Py, Pz) and its molecule ID. Out of all the species, 12 are studied, grouped into six pairs of a normal molecule and its variant:

| Normal (ID) | Variant (ID) |
|---|---|
| Carbon monoxide (211) | Carbon-13 monoxide (-211) |
| Nitric oxide (321) | Ionised NO (-321) |
| Water (2212) | Heavy water D2O (-2212) |
| Methane CH4 (3122) | Methyl ion CH3- (-3122) |
| Ethylene C2H4 (3312) | Ionised ethylene C2H3- (-3312) |
| Ozone O3 (3334) | Superoxide anion O2- (-3334) |

**Scientific question:** Are the two molecules of each pair produced at the same rate, or is one of them favoured, and does this depend on the momentum of the molecule? If there is no preference, both molecules should have the same average count within the uncertainty and the asymmetry should be zero at every momentum.

To decide whether a result is real or a statistical fluctuation, its distance from zero is expressed in units of its uncertainty (σ): below 3σ is consistent with a fluctuation, 3σ to 5σ is a hint, and 5σ or more is significant.

**Main result:** Three of the six pairs (carbon monoxide, water and methane) show a significant asymmetry (> 5σ), always in favour of the normal molecule. The other three pairs are consistent with equal production. No significant dependence of the asymmetry on momentum is observed for any pair.

## Method

**Sub-samples.** The data are split into 10 sub-samples (one per file). Every quantity is computed separately in each sub-sample. The final result is the weighted average over the sub-samples, and the **statistical uncertainty is the standard deviation (ddof = 1) of the 10 sub-sample results**, i.e. it is taken from the spread of the results.

**Empty events.** Events with 0 molecules are not counted in the number of events. The number of non-empty events per sub-sample is:

```
[461701 461835 461661 461686 461847 462025 461870 461891 461827 461650]
```

for a total of 4,617,993 events (about 7.6 % of events were empty). The largest and smallest sub-samples differ by only 0.1 % (max/min = 1.001), so the sub-samples are equivalent and their spread is a fair estimate of the statistical uncertainty.

**Quantities.**

- **Average count per event:** number of molecules of one ID / number of events, weighted over sub-samples by their number of events.
- **Asymmetry:** A = (N₁ − N₂) / (N₁ + N₂), with N₁ the count of the normal molecule and N₂ of its variant. A > 0 means the normal molecule is produced more. The overall asymmetry is the average of the sub-sample asymmetries, weighted by N₁ + N₂.
- **Significance:** |A| / uncertainty, in units of σ.
- **Momentum:** P = √(Px² + Py² + Pz²), in units of 10⁻²³ kg m/s.

**Momentum bins.** Variable-width bins are used, narrow where most molecules are and wider in the tail: 0–0.5, 0.5–1, 1–1.5, 1.5–2, 2–2.5, 2.5–3, 3–4, 4–6, 6–8 and > 8 (the last bin is open and contains all molecules above 8). Ethylene and ozone are much rarer, so adjacent bins are merged for these two pairs to keep enough molecules per bin in every sub-sample: 0–1, 1–2, 2–3, 3–6 and > 6. In every bin the asymmetry is computed per sub-sample, averaged (weighted by N₁ + N₂) and given the standard deviation of the sub-samples as uncertainty.

**Tests for momentum dependence.** For each pair the binned asymmetries are compared with:
- **zero asymmetry:** χ² = Σ (Aᵢ / σᵢ)², with the number of bins as degrees of freedom;
- **a constant equal to the overall asymmetry:** χ² = Σ ((Aᵢ − A)/σᵢ)²; if this χ² is about equal to the number of bins, the asymmetry does not depend on P;
- **a straight line:** the slope dA/dP and its uncertainty from a weighted linear fit (open tail bin excluded).

Only statistical uncertainties are evaluated. Systematic uncertainties were not studied.

## 1. What are the average counts of each molecular species and their statistical uncertainties?

| ID | Molecule | Average count per event ± stat. |
|---|---|---|
| 211 | Carbon monoxide | 19.94951 ± 0.03275 |
| -211 | Carbon-13 monoxide | 19.91721 ± 0.03187 |
| 321 | Nitric oxide | 2.50915 ± 0.00477 |
| -321 | Ionised NO | 2.50346 ± 0.00550 |
| 2212 | Water | 1.20803 ± 0.00190 |
| -2212 | Heavy water (D2O) | 1.18416 ± 0.00241 |
| 3122 | Methane (CH4) | 0.27660 ± 0.00107 |
| -3122 | Methyl ion (CH3-) | 0.27170 ± 0.00099 |
| 3312 | Ethylene (C2H4) | 0.03944 ± 0.00028 |
| -3312 | Ionised ethylene (C2H3-) | 0.03900 ± 0.00040 |
| 3334 | Ozone (O3) | 0.00119 ± 0.00004 |
| -3334 | Superoxide anion (O2-) | 0.00115 ± 0.00005 |

Carbon monoxide is by far the most common species (about 20 per event) and ozone the rarest (about 1 per 850 events). The relative uncertainty grows as the species becomes rarer: 0.2 % for carbon monoxide, 4 % for superoxide anion.

## 2. Is there any asymmetry between the normal and the variant molecule?

| Pair | Overall asymmetry ± stat. | Significance | Verdict |
|---|---|---|---|
| Carbon monoxide vs Carbon-13 monoxide | +0.00081 ± 0.00011 | 7.2 σ | significant |
| Nitric oxide vs Ionised NO | +0.00114 ± 0.00066 | 1.7 σ | not significant |
| Water vs Heavy water (D2O) | +0.00998 ± 0.00100 | 10.0 σ | significant |
| Methane (CH4) vs Methyl ion (CH3-) | +0.00894 ± 0.00106 | 8.4 σ | significant |
| Ethylene (C2H4) vs Ionised ethylene (C2H3-) | +0.0056 ± 0.0062 | 0.9 σ | not significant |
| Ozone (O3) vs Superoxide anion (O2-) | +0.015 ± 0.030 | 0.5 σ | not significant |

All six asymmetries are positive, so in every pair the normal molecule is (slightly) more frequent. Only three are significant: water (1.0 %, 10σ), methane (0.9 %, 8.4σ) and carbon monoxide (0.08 %, 7.2σ). The size of an asymmetry does not decide its significance: ozone has the largest asymmetry (1.5 %) but it is only 0.5σ from zero because ozone is so rare, while carbon monoxide has the smallest asymmetry (0.08 %) but is significant because it is so abundant that the uncertainty is tiny.

## 3. Is there any asymmetry as a function of their momentum?

### Momentum spectra

Figure 1

*Figure 1. Momentum spectrum dN/dP (counts per event per unit P, log scale) of the normal molecule (blue) and its variant (red) for each pair. Vertical error bars are the standard deviation of the 10 sub-samples, horizontal bars the bin width. The open tail bin is not drawn because it has no defined width.*

For every pair the normal and variant spectra lie on top of each other: they have the same shape, peaking at low momentum and falling steeply towards high P. Carbon monoxide and nitric oxide peak at P ≈ 0.5–1, water and methane at P ≈ 1–1.5, and ethylene and ozone at P ≈ 1–2, so heavier species are produced with somewhat higher momentum. The difference between normal and variant (at the 1 % level or below) is too small to be seen on a log scale, which is why the asymmetry is plotted directly in Figure 2.

### Asymmetry vs momentum

Figure 2

*Figure 2. Asymmetry (N₁ − N₂)/(N₁ + N₂) as a function of momentum P for each pair. Blue points: weighted mean of the sub-sample asymmetries in each bin; vertical bars: standard deviation of the 10 sub-samples (statistical uncertainty); horizontal bars: bin width. The last bin is the open tail (P > 8, or P > 6 for ethylene and ozone), drawn as wide as the previous bin. The orange band is the overall asymmetry ± its uncertainty from Question 2, the dashed line is zero asymmetry. The χ² of the points with respect to zero is printed in each panel (χ² / number of bins).*

| Pair | χ² vs zero / bins | Behaviour with momentum |
|---|---|---|
| Carbon monoxide vs Carbon-13 monoxide | 9.2 / 10 | Small and positive at all P, on the overall band; consistent with constant |
| Nitric oxide vs Ionised NO | 3.8 / 10 | Scatters around zero at all P |
| Water vs Heavy water (D2O) | 47.2 / 10 | ≈ +0.01 in every bin, flat across the whole P range |
| Methane (CH4) vs Methyl ion (CH3-) | 18.5 / 10 | ≈ +0.007 to +0.013 in every bin, flat across P |
| Ethylene (C2H4) vs Ionised ethylene (C2H3-) | 0.8 / 5 | Consistent with zero at all P (large uncertainties) |
| Ozone (O3) vs Superoxide anion (O2-) | 0.7 / 5 | Consistent with zero at all P (large uncertainties) |

The full per-bin values, the χ² against a constant asymmetry and the linear-slope fits are written by the script to `results_Q3.md`.

**Interpretation.**

- **Water and methane:** the binned asymmetries are clearly above zero (χ² vs 0 = 47.2/10 and 18.5/10) and every point lies on, or within its uncertainty of, the orange overall-asymmetry band. The asymmetry is therefore present at all momenta and has the same size everywhere: no momentum dependence is seen.
- **Carbon monoxide:** every bin is positive and follows the overall band, with the low-P bins slightly below and the 4–6 bin slightly above it, but all within their uncertainties. In a single bin the asymmetry is not significant on its own (χ² vs 0 = 9.2/10) because the statistics are divided over 10 bins, but the bins together agree with the constant 0.08 % asymmetry found in Question 2.
- **Nitric oxide, ethylene and ozone:** all bins are consistent with zero (χ² vs 0 below the number of bins), so there is no asymmetry at any momentum, in agreement with Question 2. For ethylene and ozone the uncertainties are large (up to ±0.1 for ozone at low P) because these molecules are rare, so only a large momentum-dependent asymmetry could have been detected.
- The uncertainties grow at high P (and in the first bin) for every pair because fewer molecules are produced there (see Figure 1).

**Answer:** No. Where an asymmetry exists (carbon monoxide, water, methane) it is the same at all momenta within the statistical uncertainties, and where no overall asymmetry exists (nitric oxide, ethylene, ozone) none appears in any momentum bin either.

## Conclusion

**Do we see any asymmetry? Yes, for three of the six pairs.**

- **Significant asymmetry (≥ 5σ):** water vs heavy water (+1.0 %, 10.0σ), methane vs methyl ion (+0.9 %, 8.4σ) and carbon monoxide vs carbon-13 monoxide (+0.08 %, 7.2σ). In all three cases the normal molecule is produced slightly but significantly more than its variant.
- **No significant asymmetry (< 3σ):** nitric oxide vs ionised NO (1.7σ), ethylene vs ionised ethylene (0.9σ) and ozone vs superoxide anion (0.5σ). Their production rates are equal within the statistical fluctuations of the simulation; no claim of a difference can be made.
- **Momentum dependence:** none. The momentum spectra of normal and variant molecules have the same shape, and the asymmetry of each pair is constant across the momentum range within the uncertainties.

**Limitations.** Only statistical uncertainties were evaluated; systematic effects (for example the choice of binning or the treatment of empty events) were not studied. The uncertainty used is the standard deviation of the 10 sub-samples, which describes the spread of one sub-sample; this is a conservative choice, since the uncertainty on the mean of all 10 sub-samples would be smaller by a factor √10. With only 10 sub-samples the uncertainty estimate is itself uncertain (by roughly 25 %), and for the rare species (ethylene and ozone) the analysis is limited by the small number of molecules.
