"""
rsm.py
======
Quadratic response-surface (RSM) surrogate model for AMR cooling capacity,
fit by ordinary least squares over samples spanning the design space.

The surrogate predicts

    Qc(mu0H_max, frequency, fluid_mdot, regenerator_effectiveness, span)

using linear, quadratic, and pairwise interaction terms.

Cooling capacity (Qc) is used as the response variable because it exhibits
substantial variation across the design space, making it well suited to
surrogate modeling. In contrast, the electrical COP predicted by the
current system model is comparatively insensitive to these design variables,
so constructing an RSM for COP would provide little additional value.

The fitted surrogate provides a computationally inexpensive approximation
to repeated AMRSystem.run() evaluations and is intended for design-space
exploration and multi-objective optimization.
"""

import numpy as np
from itertools import combinations
from core.mce_material import GADOLINIUM
from core.amr_cycle import AMRSystem

T_COLD_K = 291.0

VAR_NAMES = ["mu0H_max_T", "frequency_Hz", "fluid_mdot_kgs", "regen_eff", "span_K"]
BOUNDS = [[1.0, 3.0], [0.5, 5.0], [0.02, 0.2], [0.6, 0.95], [3.0, 18.0]]


def sample_design(n=300, seed=7):
    rng = np.random.default_rng(seed)
    X = np.zeros((n, 5))
    for i, (lo, hi) in enumerate(BOUNDS):
        X[:, i] = rng.uniform(lo, hi, n)
    return X


def evaluate(X):
    Y = np.zeros(len(X))
    for i, row in enumerate(X):
        mu0H, freq, mdot, eps, span = row
        sys_ = AMRSystem(material=GADOLINIUM, mu0H_max=mu0H, mass_regenerator=5.0,
                          frequency=freq, fluid_mdot=mdot, regenerator_effectiveness=eps)
        Y[i] = sys_.run(T_COLD_K, span).Qc
    return Y


def build_quadratic_features(X, names):
    n, k = X.shape
    feat_names = ["intercept"] + list(names)
    cols = [np.ones(n)] + [X[:, i] for i in range(k)]
    for i in range(k):  # squared terms
        cols.append(X[:, i] ** 2)
        feat_names.append(f"{names[i]}^2")
    for i, j in combinations(range(k), 2):  # pairwise interactions
        cols.append(X[:, i] * X[:, j])
        feat_names.append(f"{names[i]}*{names[j]}")
    return np.column_stack(cols), feat_names


def fit_rsm(out_path="results/rsm_coefficients.txt", n_train=300, n_test=100):
    X_train = sample_design(n_train, seed=7)
    Y_train = evaluate(X_train)
    X_test = sample_design(n_test, seed=99)
    Y_test = evaluate(X_test)

    Phi_train, feat_names = build_quadratic_features(X_train, VAR_NAMES)
    coeffs, *_ = np.linalg.lstsq(Phi_train, Y_train, rcond=None)

    Phi_test, _ = build_quadratic_features(X_test, VAR_NAMES)
    Y_pred = Phi_test @ coeffs
    ss_res = np.sum((Y_test - Y_pred) ** 2)
    ss_tot = np.sum((Y_test - np.mean(Y_test)) ** 2)
    r2 = 1 - ss_res / ss_tot
    rmse = np.sqrt(np.mean((Y_test - Y_pred) ** 2))

    lines = [f"Quadratic RSM surrogate for AMR cooling capacity Qc (W)",
             f"Design variables: {VAR_NAMES}",
             f"Bounds: {BOUNDS}",
             f"Train samples: {n_train}, Held-out test samples: {n_test}",
             f"Held-out R^2 = {r2:.4f}, RMSE = {rmse:.2f} W", ""]
    for name, c in zip(feat_names, coeffs):
        lines.append(f"{name:<28} {c: .6f}")

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nWrote {out_path}")
    return coeffs, feat_names, r2, rmse


if __name__ == "__main__":
    fit_rsm()
