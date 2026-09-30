"""Grandeurs conservées : énergie, quantité de mouvement, moment cinétique, barycentre.

Un bon intégrateur les conserve (Verlet : l'énergie oscille sans dériver ; Euler
explicite : elle croît sans limite — voir II.B.3 du sujet).
"""

import numpy as np

from .constantes import G as G_UNIVERSEL


def energie_cinetique(m, vit):
    return 0.5 * float(np.sum(np.asarray(m) * np.einsum("ji,ji->j", vit, vit)))


def energie_potentielle(m, pos, G=G_UNIVERSEL, adoucissement=0.0):
    """−G·Σ_{j<k} m_j·m_k / r_jk."""
    m = np.asarray(m, dtype=float)
    pos = np.asarray(pos, dtype=float)
    d = pos[None, :, :] - pos[:, None, :]
    r = np.sqrt(np.einsum("jki,jki->jk", d, d) + adoucissement ** 2)
    iu = np.triu_indices(len(m), 1)
    return float(-G * np.sum(m[iu[0]] * m[iu[1]] / r[iu]))


def energie(systeme, adoucissement=0.0):
    return (energie_cinetique(systeme.masses, systeme.vitesses)
            + energie_potentielle(systeme.masses, systeme.positions, systeme.G, adoucissement))


def quantite_de_mouvement(m, vit):
    return np.asarray(m) @ np.asarray(vit)


def moment_cinetique(m, pos, vit):
    """Σ m_j · (p_j × v_j)."""
    return np.asarray(m) @ np.cross(pos, vit)


def energies(traj, adoucissement=0.0):
    """Énergie totale à chaque instant enregistré (vecteur de longueur T)."""
    m, G = traj.systeme.masses, traj.systeme.G
    return np.array([energie_cinetique(m, v) + energie_potentielle(m, p, G, adoucissement)
                     for p, v in zip(traj.positions, traj.vitesses)])


def derive_relative(valeurs):
    """Écart relatif à la première valeur ; écart absolu si celle-ci est nulle.

    Pour une suite de vecteurs, c'est ‖x_t − x_0‖ / ‖x_0‖.
    """
    valeurs = np.asarray(valeurs, dtype=float)
    ref = valeurs[0]
    if valeurs.ndim == 1:
        return (valeurs - ref) / (abs(ref) or 1.0)
    return np.linalg.norm(valeurs - ref, axis=1) / (np.linalg.norm(ref) or 1.0)


def bilan(traj, adoucissement=0.0):
    """Dérives maximales des grandeurs conservées sur toute la trajectoire."""
    m = traj.systeme.masses
    E = energies(traj, adoucissement)
    L = np.array([moment_cinetique(m, p, v) for p, v in zip(traj.positions, traj.vitesses)])
    P = np.array([quantite_de_mouvement(m, v) for v in traj.vitesses])
    echelle_p = np.sum(m * np.linalg.norm(traj.vitesses[0], axis=1)) or 1.0
    return {
        "energie": float(np.max(np.abs(derive_relative(E)))),
        "moment_cinetique": float(np.max(derive_relative(L))),
        "quantite_de_mouvement": float(np.max(np.linalg.norm(P - P[0], axis=1)) / echelle_p),
    }


def periode_moyenne(t, x):
    """Période d'un signal oscillant, estimée par ses passages ascendants par la moyenne.

    Retourne None si moins de deux passages sont observés.
    """
    x = np.asarray(x, dtype=float) - np.mean(x)
    idx = np.where((x[:-1] < 0) & (x[1:] >= 0))[0]
    if len(idx) < 2:
        return None
    # interpolation linéaire du passage par zéro
    tz = t[idx] - x[idx] * (t[idx + 1] - t[idx]) / (x[idx + 1] - x[idx])
    return float(np.mean(np.diff(tz)))
