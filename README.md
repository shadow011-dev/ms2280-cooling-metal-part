# MS2280 Project 7: Cooling a metal part

## Files
- `Report.pdf` - the report (definitions, assumptions, derivation, results, checks, limitations, references).
- `cooling-lab.html` - interactive web model. Open it in a recent Chrome, Edge, Firefox or Safari (internet needed for fonts, Three.js and Chart.js). Every graph on the page comes from one model; the Measurements section fits a time constant to pasted time/temperature data.
- `cooling_model.py` - Python model that reproduces every number and figure in the report.
- `results.csv` - the nine cases (3 sizes x 3 cooling conditions): tau, Biot number, centre-surface gap, verdict.
- `figures/` - plots written by the script.

## Reproduce the results
1. Install Python 3 with `numpy` and `matplotlib`.
2. Run `python3 cooling_model.py` in this folder.
3. It prints the nine-case table, the energy balance, the metal comparison and the limiting cases, and rewrites `results.csv` and `figures/`.

## Model summary
Lumped capacitance: m c dT/dt = -h A (T - Tinf), so tau = rho c V / (h A) and T = Tinf + (T0 - Tinf) exp(-t/tau). Valid when Bi = h Lc / k < 0.1, with Lc = V/A = w/4 for a long square bar. The spatial check solves 1D transient conduction by implicit finite differences in a sphere with the same V/A (R = 0.75 w).
Test part: steel bar, 300 mm long, 12.5 / 25 / 50 mm square, T0 = 900 C, Tinf = 25 C; h = 10 (still air), 100 (fan), 1000 W/m2K (water).
