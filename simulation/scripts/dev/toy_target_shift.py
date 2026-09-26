"""Exhaustive three-event example: the population minimiser of the DEBM consensus objective
changes with the stage composition even with ONE common measurement model and identical
conditional laws in every group.

Model. Events A, B, C occur in the true order A < B < C. A participant at stage k has the first k events
occurred. Each biomarker is N(0,1) before and N(delta_e,1) after its event; the (common, pooled)
measurement model turns a value x into q_e(x) = pi phi(x-delta_e) / (pi phi(x-delta_e) + (1-pi) phi(x)),
with a common mixing weight pi. The participant-level loss is pyebm's weighted Kendall discordance
   l(S; q) = sum over pairs (a before b in S, q_a < q_b) of (q_b - q_a),
and the best-fitting ordering at composition w = (w_0,...,w_3) over stages is argmin_S sum_k w_k L_k(S),
L_k(S) = E[l(S;q) | stage k]. L_k is computed by Monte Carlo with a common random stream per stage.
"""
import itertools, numpy as np
from math import erf, sqrt
rng = np.random.default_rng(20260914)
EV = ['A', 'B', 'C']; ORDERS = list(itertools.permutations(range(3)))
def phi(x): return np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)
def auc(d): return 0.5 * (1 + erf(d / 2))
def posterior(x, delta, pi): return pi * phi(x - delta) / (pi * phi(x - delta) + (1 - pi) * phi(x))
def loss(S, q):
    """weighted Kendall discordance between ordering S (tuple of event indices) and posterior vector q (n x 3)"""
    tot = np.zeros(q.shape[0])
    for i in range(3):
        for j in range(i + 1, 3):
            a, b = S[i], S[j]                     # S says a before b
            d = q[:, b] - q[:, a]; tot += np.where(d > 0, d, 0.)
    return tot
def stage_losses(delta, pi, n=400_000):
    L = np.zeros((4, len(ORDERS)))
    for k in range(4):
        occurred = np.arange(3) < k
        x = rng.standard_normal((n, 3)) + occurred * delta
        q = np.column_stack([posterior(x[:, e], delta[e], pi) for e in range(3)])
        for s, S in enumerate(ORDERS): L[k, s] = loss(S, q).mean()
    return L
def name(S): return '<'.join(EV[e] for e in S)
for delta, label in (((1.0, 2.5, 1.5), 'A poorly separable (AUC .76), B sharp (.96), C middling (.86)'),
                     ((1.5, 1.5, 1.5), 'equal separability (AUC .86)')):
    delta = np.array(delta); pi = 0.5
    L = stage_losses(delta, pi)
    print(f"\n=== {label}; common measurement model, mixing weight {pi}")
    print("per-stage minimiser:", {k: name(ORDERS[int(np.argmin(L[k]))]) for k in range(4)})
    print("stage 0 losses:", {name(S): round(L[0, s], 4) for s, S in enumerate(ORDERS)})
    print("stage 3 losses:", {name(S): round(L[3, s], 4) for s, S in enumerate(ORDERS)})
    # compositions: early-heavy vs late-heavy, sweep the late fraction
    print("composition sweep (w over stages 0..3):")
    last = None
    for f in np.linspace(0, 1, 21):
        w = np.array([(1 - f) * 0.6, (1 - f) * 0.4, f * 0.4, f * 0.6])   # early mass on stages 0/1, late mass on 2/3
        obj = w @ L; S = ORDERS[int(np.argmin(obj))]
        if name(S) != last: print(f"  late fraction {f:.2f}: target {name(S)}   (objective {obj.min():.4f}; truth A<B<C costs {obj[0]:.4f})"); last = name(S)
    # the two ADNI-like compositions from Table 2 mapped to stages: CN -> stages 0/1, MCI -> 1/2, AD -> 2/3
    for grp, (cn, mci, ad) in (('e2', (57, 6, 12)), ('e4', (110, 156, 219))):
        tot = cn + mci + ad; w = np.array([cn * .7, cn * .3 + mci * .5, mci * .5 + ad * .3, ad * .7]) / tot
        obj = w @ L; print(f"  {grp} composition: target {name(ORDERS[int(np.argmin(obj))])}  losses {dict((name(S), round(o,4)) for S, o in zip(ORDERS, obj))}")
