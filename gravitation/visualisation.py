"""Graphiques (matplotlib). Même palette que le projet algoTrade."""

import numpy as np

from .diagnostics import derive_relative, energies

BLEU, ORANGE, VERT, ROUGE, GRIS, ENCRE = ("#2a78d6", "#eb6834", "#008300", "#e34948",
                                          "#8a8983", "#0b0b0b")
COULEURS = [ORANGE, BLEU, VERT, ROUGE, "#9b59b6", "#8c6d31", "#17a2b8", GRIS, ENCRE]


def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#eeeeea")


def tracer_orbites(traj, unite=1.0, etiquette="m", titre="", limite=None):
    """Trajectoires dans le plan xy ; `limite` restreint la fenêtre (en `unite`)."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 7))
    for j, nom in enumerate(traj.systeme.noms):
        xy = traj.positions[:, j, :2] / unite
        c = COULEURS[j % len(COULEURS)]
        ax.plot(xy[:, 0], xy[:, 1], color=c, lw=1.2, label=nom)
        ax.scatter(*xy[-1], color=c, s=25, zorder=3)
    ax.set_aspect("equal")
    if limite:
        ax.set_xlim(-limite, limite)
        ax.set_ylim(-limite, limite)
    ax.set_xlabel(f"x ({etiquette})")
    ax.set_ylabel(f"y ({etiquette})")
    ax.legend(loc="upper right", frameon=False, fontsize=8, ncol=2)
    ax.set_title(titre or "Trajectoires", fontweight="bold")
    _style(ax)
    fig.tight_layout()
    return fig


def tracer_energie(trajectoires, unite_temps=1.0, etiquette="s", titre=""):
    """Dérive relative de l'énergie totale pour plusieurs schémas (échelle log symétrique)."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for k, traj in enumerate(trajectoires):
        d = derive_relative(energies(traj))
        ax.plot(traj.t / unite_temps, d, color=COULEURS[(k + 1) % len(COULEURS)], lw=1.6,
                label=traj.methode)
    ax.set_yscale("symlog", linthresh=1e-12)
    ax.axhline(0, color="#d9d8d3", lw=1)
    ax.set_xlabel(f"temps ({etiquette})")
    ax.set_ylabel("(E − E₀) / |E₀|")
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(titre or "Conservation de l'énergie", fontweight="bold")
    _style(ax)
    fig.tight_layout()
    return fig


def tracer_oscillateur(t, euler_yz, verlet_yz, exact_y, omega):
    """Portrait de phase (y, z) et énergie : Euler diverge, Verlet reste sur l'ellipse."""
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    th = np.linspace(0, 2 * np.pi, 400)
    y0 = exact_y[0]
    ax1.plot(y0 * np.cos(th), -y0 * omega * np.sin(th), color="#d9d8d3", lw=6, label="exact")
    for (nom, (y, z)), c in zip((("Euler", euler_yz), ("Verlet", verlet_yz)), (ROUGE, BLEU)):
        ax1.plot(y, z, color=c, lw=1.3, label=nom)
        e = 0.5 * np.array(z) ** 2 + 0.5 * omega ** 2 * np.array(y) ** 2
        ax2.plot(t, e / e[0], color=c, lw=1.8, label=nom)
    ax1.set_xlabel("y")
    ax1.set_ylabel("z = y'")
    ax1.legend(frameon=False)
    ax1.set_title("Portrait de phase", fontweight="bold")
    ax2.set_xlabel("t")
    ax2.set_ylabel("E / E₀")
    ax2.legend(frameon=False)
    ax2.set_title("Énergie", fontweight="bold")
    for ax in (ax1, ax2):
        _style(ax)
    fig.tight_layout()
    return fig


def tracer_complexite(tailles, mesures, pentes):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for (nom, temps), c in zip(mesures.items(), (ROUGE, BLEU)):
        ax.loglog(tailles, temps, "o-", color=c, lw=1.8, label=f"{nom} (pente {pentes[nom]:.2f})")
    ax.set_xlabel("N (nombre de corps)")
    ax.set_ylabel("temps d'un pas (s)")
    ax.legend(frameon=False)
    ax.set_title("Complexité de etat_suiv", fontweight="bold")
    _style(ax)
    fig.tight_layout()
    return fig
