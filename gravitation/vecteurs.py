"""Opérations vectorielles sur des listes (I.A du sujet) : smul, vsom, vdif."""

import math


def smul(scalaire, lst):
    """Multiplie chaque élément de `lst` par `scalaire`."""
    return [scalaire * x for x in lst]


def vsom(lst1, lst2):
    """Somme terme à terme de deux listes de même longueur."""
    _verifier(lst1, lst2)
    return [x + y for x, y in zip(lst1, lst2)]


def vdif(lst1, lst2):
    """Différence terme à terme `lst1 - lst2`."""
    _verifier(lst1, lst2)
    return [x - y for x, y in zip(lst1, lst2)]


def produit_scalaire(lst1, lst2):
    _verifier(lst1, lst2)
    return sum(x * y for x, y in zip(lst1, lst2))


def norme(lst):
    return math.sqrt(sum(x * x for x in lst))


def _verifier(lst1, lst2):
    if len(lst1) != len(lst2):
        raise ValueError(f"listes de longueurs différentes ({len(lst1)} et {len(lst2)})")
