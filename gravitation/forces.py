"""Forces gravitationnelles (partie III.A du sujet).

Deux interfaces au même calcul :

* `force2` / `forceN` : listes Python, fidèles au sujet, faciles à lire ;
* `accelerations` : version numpy vectorisée qui calcule d'un coup les N accélérations.
"""

import math

import numpy as np

from .constantes import G as G_UNIVERSEL


def force2(m1, p1, m2, p2, G=G_UNIVERSEL):
    """Force exercée par le corps 2 sur le corps 1 : G·m1·m2 / r³ · (p2 − p1)."""
    r = [p2[i] - p1[i] for i in range(3)]
    d = math.sqrt(sum(c * c for c in r))
    if d == 0.0:
        raise ValueError("deux corps à la même position : force infinie")
    coeff = G * m1 * m2 / d ** 3
    return [coeff * c for c in r]


def forceN(j, m, pos, G=G_UNIVERSEL):
    """Force totale exercée sur le corps j par tous les autres."""
    total = [0.0, 0.0, 0.0]
    for k in range(len(m)):
        if k != j:
            f = force2(m[j], pos[j], m[k], pos[k], G)
            total = [total[i] + f[i] for i in range(3)]
    return total


def accelerations(m, pos, G=G_UNIVERSEL, adoucissement=0.0):
    """Accélérations des N corps, sans boucle Python : a_j = G·Σ_k m_k·(p_k − p_j)/r_jk³.

    `adoucissement` (m) ajoute ε² à r² : évite les forces infinies lors des
    rencontres très proches (0 par défaut = gravitation exacte).
    """
    m = np.asarray(m, dtype=float)
    pos = np.asarray(pos, dtype=float)
    d = pos[None, :, :] - pos[:, None, :]            # d[j, k] = p_k − p_j
    r2 = np.einsum("jki,jki->jk", d, d) + adoucissement ** 2
    np.fill_diagonal(r2, 1.0)                        # évite 0⁻¹·⁵ sur la diagonale
    if np.any(r2 == 0.0):
        raise ValueError("deux corps à la même position : force infinie")
    poids = m[None, :] * r2 ** -1.5
    np.fill_diagonal(poids, 0.0)                     # pas d'auto-attraction
    return G * np.einsum("jk,jki->ji", poids, d)
