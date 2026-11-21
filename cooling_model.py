"""MS2280 - Cooling a metal part: lumped capacitance model with independent checks.

Run:  python3 cooling_model.py
Needs: numpy, matplotlib.  Writes results.csv and figures/*.png next to this file.
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
METALS = {"Steel": (7850, 490, 45), "Copper": (8960, 385, 401), "Aluminium": (2700, 900, 237)}  # rho, c, k
RHO, C, K = METALS["Steel"]
T0, TINF, LBAR = 900.0, 25.0, 0.3
TH0 = T0 - TINF
SIZES = [0.0125, 0.025, 0.05]                               # bar side w (m)
CONDS = [("Still air", 10), ("Fan", 100), ("Water", 1000)]  # h (W/m^2 K)


def tau(w, h, rho=RHO, c=C):          # tau = rho c V / (h A) = rho c Lc / h, Lc = w/4
    return rho * c * (w / 4) / h


def biot(w, h, k=K):                  # Bi = h Lc / k
    return h * (w / 4) / k


def T_lumped(t, w, h):                # T = Tinf + theta0 exp(-t/tau)
    return TINF + TH0 * np.exp(-t / tau(w, h))


def spatial(w, h, n=240, tmax_tau=4, N=24):
    """1D transient conduction, implicit finite differences, in a sphere with the
    same V/A as the bar section (R = 0.75 w) and a convective surface.
    Returns t, centre excess temperature, surface excess temperature."""
    R = 0.75 * w
    dr, tmax = R / N, tmax_tau * tau(w, h)
    f = K / (RHO * C) * (tmax / n) / dr**2
    b = h / K
    lo, di, up = np.zeros(N + 1), np.zeros(N + 1), np.zeros(N + 1)
    for i in range(N + 1):
        if i == 0:
            di[i], up[i] = 1 + 6 * f, -6 * f
        elif i == N:
            lo[i], di[i] = -2 * f, 1 + f * (2 + 2 * dr * b + 2 * b * dr / N)
        else:
            lo[i], di[i], up[i] = -f * (1 - 1 / i), 1 + 2 * f, -f * (1 + 1 / i)
    u = np.full(N + 1, TH0)
    t, ctr, srf = [0.0], [TH0], [TH0]
    for s in range(n):                # Thomas algorithm for the tridiagonal system
        cp, dp = np.zeros(N + 1), np.zeros(N + 1)
        cp[0], dp[0] = up[0] / di[0], u[0] / di[0]
        for i in range(1, N + 1):
            d = di[i] - lo[i] * cp[i - 1]
            cp[i], dp[i] = up[i] / d, (u[i] - lo[i] * dp[i - 1]) / d
        u[N] = dp[N]
        for i in range(N - 1, -1, -1):
            u[i] = dp[i] - cp[i] * u[i + 1]
        t.append((s + 1) * tmax / n); ctr.append(u[0]); srf.append(u[N])
    return np.array(t), np.array(ctr), np.array(srf)


def max_gap(w, h):                    # largest (centre - surface) / initial difference
    _, c, s = spatial(w, h)
    return ((c - s) / TH0).max()


def energy_check(w=0.025, h=100, n=300):
    ta, m, A = tau(w, h), RHO * w * w * LBAR, 4 * w * LBAR
    t = np.linspace(0, 6 * ta, n + 1)
    q = A * h * TH0 * np.exp(-t / ta)
    num = np.sum((q[1:] + q[:-1]) / 2 * np.diff(t))
    ana = m * C * TH0 * (1 - np.exp(-6))
    return m, A, ta, num / 1e3, ana / 1e3, abs(num - ana) / ana * 100


def main():
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    rows = []
    print(f"{'size':>8} {'condition':>10} {'h':>5} {'tau (s)':>9} {'Bi':>7} {'gap %':>6}  verdict")
    for w in SIZES:
        for name, h in CONDS:
            bi, g = biot(w, h), max_gap(w, h) * 100
            verdict = "Lumped valid" if bi < 0.1 else "Marginal"
            rows.append([w * 1e3, name, h, round(tau(w, h), 1), round(bi, 4), round(g, 1), verdict])
            print(f"{w*1e3:>6g}mm {name:>10} {h:>5} {tau(w,h):>9.1f} {bi:>7.4f} {g:>6.1f}  {verdict}")
    with open(os.path.join(HERE, "results.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["size_mm", "condition", "h_W_m2K", "tau_s", "Biot", "max_centre_surface_gap_pct", "verdict"])
        wr.writerows(rows)

    m, A, ta, num, ana, err = energy_check()
    print(f"\nEnergy balance (25 mm, h=100): m={m:.3f} kg, A={A:.3f} m2, tau={ta:.0f} s")
    print(f"  integral q dt = {num:.3f} kJ, m c theta0 (1-e^-6) = {ana:.3f} kJ, error = {err:.1e} %")
    print("\nMetals (25 mm, h=100):")
    for name, (r, c, k) in METALS.items():
        print(f"  {name:10s} tau = {tau(.025,100,r,c)/60:.1f} min, Bi = {biot(.025,100,k):.4f}")
    print(f"\nLimits: T(0) = {T_lumped(0,.025,100):.1f} C; 6 tau -> within {TH0*np.exp(-6):.2f} K of Tinf; "
          f"h=0.01 -> tau = {tau(.025,.01)/3600:.1f} h; h=1e5 -> tau = {tau(.025,1e5):.2f} s")

    # figures: temperature, energy / q / ln theta, lumped vs spatial
    w, h = 0.025, 100
    t = np.linspace(0, 40 * 60, 400)
    fig, ax = plt.subplots(figsize=(5, 3))
    for s in SIZES:
        ax.plot(t / 60, T_lumped(t, s, h), label=f"{s*1e3:g} mm, tau = {tau(s,h)/60:.1f} min")
    ax.set(xlabel="time (min)", ylabel="temperature (C)", title="Three sizes, h = 100"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "sizes.png"), dpi=200); plt.close()

    t = np.linspace(0, 6 * ta, 301)
    E = m * C * TH0 * np.exp(-t / ta) / 1e3
    q = A * h * TH0 * np.exp(-t / ta)
    fig, ax = plt.subplots(1, 3, figsize=(10, 3))
    ax[0].plot(t / 60, E); ax[0].set(xlabel="time (min)", ylabel="E (kJ)", title="Stored energy")
    ax[1].plot(t / 60, q); ax[1].set(xlabel="time (min)", ylabel="q (W)", title="Heat transfer rate")
    ax[2].plot(t, np.log(TH0 * np.exp(-t / ta))); ax[2].set(xlabel="time (s)", ylabel="ln theta", title="ln theta vs t (slope -1/tau)")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "quantities.png"), dpi=200); plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(9, 3), sharey=True)
    for a, (s, hh) in zip(ax, [(0.025, 10), (0.05, 1000)]):
        tt, c, sf = spatial(s, hh); x = tt / tau(s, hh)
        a.plot(x, -x, "--k", label="lumped"); a.plot(x, np.log(c / TH0), label="centre"); a.plot(x, np.log(sf / TH0), label="surface")
        a.set(xlabel="t / tau", ylim=(-4, 0.1), title=f"{s*1e3:g} mm, h = {hh}, Bi = {biot(s,hh):.3f}"); a.legend()
    ax[0].set_ylabel("ln theta")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "spatial_check.png"), dpi=200); plt.close()
    print("\nWrote results.csv and figures/")


if __name__ == "__main__":
    main()
