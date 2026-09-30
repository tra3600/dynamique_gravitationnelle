"""Équation différentielle y''(t) = f(y(t)) : schémas d'Euler et de Verlet (partie II du sujet).

On passe au système du premier ordre y' = z, z' = f(y). `y` peut être un flottant
ou un tableau numpy (le même code sert à un oscillateur ou à un corps en 3D).
"""

import math


def grille(tmin, tmax, n):
    """Les n instants t_i = tmin + i·h, avec h = (tmax − tmin)/(n − 1)."""
    _verifier(tmin, tmax, n)
    h = (tmax - tmin) / (n - 1)
    return [tmin + i * h for i in range(n)]


def euler(f, y0, z0, tmin, tmax, n):
    """Euler explicite : y_{i+1} = y_i + h·z_i,  z_{i+1} = z_i + h·f(y_i).

    Retourne les deux listes (y_i) et (z_i), de longueur n.
    """
    _verifier(tmin, tmax, n)
    h = (tmax - tmin) / (n - 1)
    y, z = [y0], [z0]
    for _ in range(n - 1):
        y.append(y[-1] + h * z[-1])
        z.append(z[-1] + h * f(y[-2]))
    return y, z


def verlet(f, y0, z0, tmin, tmax, n):
    """Verlet (vitesses) : y_{i+1} = y_i + h·z_i + h²/2·f(y_i),
    z_{i+1} = z_i + h/2·(f(y_i) + f(y_{i+1})).

    Une seule évaluation de f par pas : f(y_{i+1}) est réutilisée au pas suivant.
    """
    _verifier(tmin, tmax, n)
    h = (tmax - tmin) / (n - 1)
    y, z = [y0], [z0]
    fi = f(y0)
    for _ in range(n - 1):
        y_suiv = y[-1] + h * z[-1] + 0.5 * h * h * fi
        f_suiv = f(y_suiv)
        z.append(z[-1] + 0.5 * h * (fi + f_suiv))
        y.append(y_suiv)
        fi = f_suiv
    return y, z


def oscillateur(omega):
    """Le second membre f(y) = −ω²·y de l'oscillateur harmonique."""
    return lambda y: -omega * omega * y


def energie_oscillateur(y, z, omega):
    """E = ½·z² + ½·ω²·y², constante pour la vraie solution."""
    return [0.5 * zi * zi + 0.5 * omega * omega * yi * yi for yi, zi in zip(y, z)]


def _verifier(tmin, tmax, n):
    if n < 2:
        raise ValueError("n doit être au moins 2")
    if not tmax > tmin:
        raise ValueError("il faut tmax > tmin")
    if not math.isfinite(tmin) or not math.isfinite(tmax):
        raise ValueError("bornes non finies")
