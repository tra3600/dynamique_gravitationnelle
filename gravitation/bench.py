"""Mesure de complexité de etat_suiv (III.B.4–5 du sujet) : temps en fonction de N."""

import time

import numpy as np

from .forces import accelerations
from .integrateurs import etat_suiv
from .systemes import amas


def _meilleur_temps(fonction, repetitions):
    meilleur = float("inf")
    for _ in range(repetitions):
        debut = time.perf_counter()
        fonction()
        meilleur = min(meilleur, time.perf_counter() - debut)
    return meilleur


def mesurer(tailles, repetitions=3, graine=0, h=60.0):
    """Temps (s) d'un pas pour chaque N : version listes (sujet) et version numpy."""
    listes, vectorise = [], []
    for n in tailles:
        s = amas(n, graine)
        m, p, v = s.en_listes()
        listes.append(_meilleur_temps(lambda: etat_suiv(m, p, v, h, s.G), repetitions))
        vectorise.append(_meilleur_temps(lambda: accelerations(s.masses, s.positions, s.G),
                                         repetitions))
    return {"listes": listes, "numpy": vectorise}


def pente_loglog(tailles, temps):
    """Pente a de la droite ln(τ) = a·ln(N) + b : le temps croît comme N^a."""
    a, _ = np.polyfit(np.log(tailles), np.log(temps), 1)
    return float(a)
