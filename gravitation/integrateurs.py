"""Intégration du problème à N corps (partie III.B du sujet et au-delà).

* `pos_suiv` / `etat_suiv` : les fonctions du sujet, sur des listes Python ;
* `simuler` : boucle complète, vectorisée, avec le choix du schéma
  (`euler`, `euler-symplectique`, `verlet`, `rk4`) ;
* `simuler_reference` : solution de référence à pas adaptatif (scipy, DOP853).
"""

from dataclasses import dataclass

import numpy as np

from .constantes import G as G_UNIVERSEL
from .forces import accelerations, forceN


# ----------------------------------------------------------------------
# Fonctions du sujet (listes Python)
def pos_suiv(m, pos, vit, h, G=G_UNIVERSEL):
    """Positions à t_{i+1} : p + h·v + h²/2·F/m."""
    nouvelles = []
    for j in range(len(m)):
        F = forceN(j, m, pos, G)
        nouvelles.append([pos[j][k] + h * vit[j][k] + 0.5 * h * h * F[k] / m[j]
                          for k in range(3)])
    return nouvelles


def etat_suiv(m, pos, vit, h, G=G_UNIVERSEL):
    """(positions, vitesses) à t_{i+1} par le schéma de Verlet.

    Les forces à t_i et à t_{i+1} sont calculées une seule fois par corps.
    """
    nouvelles = pos_suiv(m, pos, vit, h, G)
    vitesses = []
    for j in range(len(m)):
        F, F_suiv = forceN(j, m, pos, G), forceN(j, m, nouvelles, G)
        vitesses.append([vit[j][k] + 0.5 * h * (F[k] + F_suiv[k]) / m[j] for k in range(3)])
    return nouvelles, vitesses


# ----------------------------------------------------------------------
# Schémas vectorisés : chacun reçoit f(pos) -> accélérations et renvoie
# (pos, vit, acc) où acc peut être None si le schéma ne la réutilise pas.
def _euler(f, p, v, a, h):
    return p + h * v, v + h * (f(p) if a is None else a), None


def _euler_symplectique(f, p, v, a, h):
    v = v + h * (f(p) if a is None else a)
    return p + h * v, v, None


def _verlet(f, p, v, a, h):
    a = f(p) if a is None else a
    p_suiv = p + h * v + 0.5 * h * h * a
    a_suiv = f(p_suiv)
    return p_suiv, v + 0.5 * h * (a + a_suiv), a_suiv


def _rk4(f, p, v, a, h):
    a1 = f(p) if a is None else a
    p2, v2 = p + 0.5 * h * v, v + 0.5 * h * a1
    a2 = f(p2)
    p3, v3 = p + 0.5 * h * v2, v + 0.5 * h * a2
    a3 = f(p3)
    p4, v4 = p + h * v3, v + h * a3
    a4 = f(p4)
    return (p + h / 6 * (v + 2 * v2 + 2 * v3 + v4),
            v + h / 6 * (a1 + 2 * a2 + 2 * a3 + a4), None)


METHODES = {"euler": _euler, "euler-symplectique": _euler_symplectique,
            "verlet": _verlet, "rk4": _rk4}


@dataclass
class Trajectoire:
    """Résultat d'une simulation : t (T,), positions (T, N, 3), vitesses (T, N, 3)."""

    t: np.ndarray
    positions: np.ndarray
    vitesses: np.ndarray
    systeme: object          # état initial (masses, noms, G)
    methode: str = ""
    duree_calcul: float = 0.0

    def etat_final(self):
        s = self.systeme.copie()
        s.positions, s.vitesses = self.positions[-1].copy(), self.vitesses[-1].copy()
        return s


def simuler(systeme, duree, h, methode="verlet", enregistrer_tous=1, adoucissement=0.0):
    """Fait évoluer `systeme` pendant `duree` secondes avec le pas `h`.

    `enregistrer_tous` : ne garde qu'un pas sur k (économise la mémoire sur de longues
    simulations). Le dernier instant est toujours conservé.
    """
    import time

    if methode not in METHODES:
        raise ValueError(f"méthode inconnue : {methode!r} (choix : {', '.join(METHODES)})")
    if h <= 0 or duree <= 0:
        raise ValueError("duree et h doivent être strictement positifs")
    pas = METHODES[methode]
    n = max(1, int(round(duree / h)))
    m, G = systeme.masses, systeme.G

    def f(p):
        return accelerations(m, p, G, adoucissement)

    p, v, a = systeme.positions.copy(), systeme.vitesses.copy(), None
    ts, ps, vs = [0.0], [p], [v]
    debut = time.perf_counter()
    for i in range(1, n + 1):
        p, v, a = pas(f, p, v, a, h)
        if i % enregistrer_tous == 0 or i == n:
            ts.append(i * h)
            ps.append(p)
            vs.append(v)
    return Trajectoire(np.array(ts), np.array(ps), np.array(vs), systeme, methode,
                       time.perf_counter() - debut)


def simuler_reference(systeme, duree, nb_points=500, rtol=1e-11, adoucissement=0.0):
    """Solution de référence : Dormand–Prince d'ordre 8 à pas adaptatif (scipy)."""
    try:
        from scipy.integrate import solve_ivp
    except ImportError as e:  # pragma: no cover
        raise ImportError("scipy est nécessaire pour la solution de référence") from e
    import time

    n, m, G = systeme.n, systeme.masses, systeme.G

    def second_membre(_, y):
        p, v = y[:3 * n].reshape(n, 3), y[3 * n:].reshape(n, 3)
        return np.concatenate([v.ravel(), accelerations(m, p, G, adoucissement).ravel()])

    y0 = np.concatenate([systeme.positions.ravel(), systeme.vitesses.ravel()])
    echelle = max(np.abs(systeme.positions).max(), 1.0)
    t = np.linspace(0.0, duree, nb_points)
    debut = time.perf_counter()
    sol = solve_ivp(second_membre, (0.0, duree), y0, method="DOP853", t_eval=t,
                    rtol=rtol, atol=rtol * echelle)
    if not sol.success:  # pragma: no cover
        raise RuntimeError(sol.message)
    y = sol.y.T
    return Trajectoire(sol.t, y[:, :3 * n].reshape(-1, n, 3), y[:, 3 * n:].reshape(-1, n, 3),
                       systeme, "référence (DOP853)", time.perf_counter() - debut)
