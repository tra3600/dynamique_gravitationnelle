"""Systèmes de N corps et scénarios prêts à l'emploi (unités SI)."""

import math
from dataclasses import dataclass, field

import numpy as np

from .constantes import G as G_UNIVERSEL, JOUR, UA


@dataclass
class Systeme:
    """N corps : masses (N,), positions (N, 3) en m, vitesses (N, 3) en m/s."""

    noms: list
    masses: np.ndarray
    positions: np.ndarray
    vitesses: np.ndarray
    G: float = G_UNIVERSEL

    def __post_init__(self):
        self.masses = np.asarray(self.masses, dtype=float)
        self.positions = np.asarray(self.positions, dtype=float)
        self.vitesses = np.asarray(self.vitesses, dtype=float)
        n = len(self.masses)
        if self.positions.shape != (n, 3) or self.vitesses.shape != (n, 3):
            raise ValueError("positions et vitesses doivent être de forme (N, 3)")
        if len(self.noms) != n:
            raise ValueError("un nom par corps est attendu")
        if np.any(self.masses <= 0):
            raise ValueError("les masses doivent être strictement positives")

    @property
    def n(self):
        return len(self.masses)

    def copie(self):
        return Systeme(list(self.noms), self.masses.copy(), self.positions.copy(),
                       self.vitesses.copy(), self.G)

    def barycentre(self):
        return self.masses @ self.positions / self.masses.sum()

    def vitesse_barycentre(self):
        return self.masses @ self.vitesses / self.masses.sum()

    def dans_le_referentiel_barycentrique(self):
        """Copie recentrée : barycentre à l'origine et quantité de mouvement nulle."""
        s = self.copie()
        s.positions -= s.barycentre()
        s.vitesses -= s.vitesse_barycentre()
        return s

    def en_listes(self):
        """(masses, positions, vitesses) en listes Python, pour les fonctions du sujet."""
        return self.masses.tolist(), self.positions.tolist(), self.vitesses.tolist()


def orbite_circulaire(M, m, a, angle=0.0, G=G_UNIVERSEL):
    """Position et vitesse relatives d'un corps de masse m sur une orbite circulaire
    de rayon a autour d'un corps de masse M (plan xy)."""
    v = math.sqrt(G * (M + m) / a)
    return (np.array([a * math.cos(angle), a * math.sin(angle), 0.0]),
            np.array([-v * math.sin(angle), v * math.cos(angle), 0.0]))


def terre_lune():
    m_t, m_l, a = 5.972e24, 7.348e22, 3.844e8
    p, v = orbite_circulaire(m_t, m_l, a)
    return Systeme(["Terre", "Lune"], [m_t, m_l], [[0, 0, 0], p], [[0, 0, 0], v]
                   ).dans_le_referentiel_barycentrique()


def soleil_terre():
    m_s, m_t = 1.98847e30, 5.97217e24
    p, v = orbite_circulaire(m_s, m_t, UA)
    return Systeme(["Soleil", "Terre"], [m_s, m_t], [[0, 0, 0], p], [[0, 0, 0], v]
                   ).dans_le_referentiel_barycentrique()


# nom, masse (kg), demi-grand axe (UA), phase initiale (rad)
_PLANETES = [
    ("Mercure", 3.3011e23, 0.387098, 0.0),
    ("Vénus", 4.8675e24, 0.723332, 0.9),
    ("Terre", 5.97217e24, 1.000000, 1.9),
    ("Mars", 6.4171e23, 1.523679, 2.8),
    ("Jupiter", 1.89813e27, 5.2044, 3.6),
    ("Saturne", 5.6834e26, 9.5826, 4.4),
    ("Uranus", 8.6810e25, 19.2184, 5.2),
    ("Neptune", 1.02413e26, 30.110, 6.0),
]


def systeme_solaire(planetes=None):
    """Le Soleil et les huit planètes sur des orbites circulaires coplanaires.

    `planetes` : noms à conserver (par défaut toutes).
    """
    m_s = 1.98847e30
    noms, masses, pos, vit = ["Soleil"], [m_s], [np.zeros(3)], [np.zeros(3)]
    for nom, m, a_ua, phase in _PLANETES:
        if planetes is not None and nom not in planetes:
            continue
        p, v = orbite_circulaire(m_s, m, a_ua * UA, phase)
        noms.append(nom)
        masses.append(m)
        pos.append(p)
        vit.append(v)
    return Systeme(noms, masses, pos, vit).dans_le_referentiel_barycentrique()


def huit():
    """Chorégraphie en « huit » de trois corps égaux (Chenciner–Montgomery), G = 1.

    Solution périodique exacte, de période ≈ 6,32591398 : un excellent test d'intégrateur.
    """
    x, y = 0.97000436, -0.24308753
    vx, vy = -0.93240737, -0.86473146
    return Systeme(["A", "B", "C"], [1.0, 1.0, 1.0],
                   [[-x, -y, 0], [x, y, 0], [0, 0, 0]],
                   [[-vx / 2, -vy / 2, 0], [-vx / 2, -vy / 2, 0], [vx, vy, 0]], G=1.0)


PERIODE_HUIT = 6.32591398


def amas(n, graine=0, G=G_UNIVERSEL):
    """N corps de masses et positions aléatoires, au repos dans le référentiel barycentrique
    (sert surtout aux mesures de complexité)."""
    rng = np.random.default_rng(graine)
    masses = rng.uniform(1e22, 1e24, n)
    pos = rng.uniform(-1, 1, (n, 3)) * UA
    return Systeme([f"c{i}" for i in range(n)], masses, pos, rng.normal(0, 1e3, (n, 3)), G
                   ).dans_le_referentiel_barycentrique()


@dataclass
class Scenario:
    """Un système avec ses réglages par défaut (durée et pas exprimés dans `unite`)."""

    nom: str
    description: str
    construire: object
    duree: float
    pas: float
    unite: float = JOUR       # secondes par unité de temps
    etiquette: str = "jours"  # nom de l'unité
    unite_longueur: float = UA
    etiquette_longueur: str = "UA"
    options: dict = field(default_factory=dict)


SCENARIOS = {
    "terre-lune": Scenario("terre-lune", "la Lune en orbite circulaire autour de la Terre",
                           terre_lune, duree=27.3, pas=0.05,
                           unite_longueur=1e8, etiquette_longueur="10⁸ m"),
    "soleil-terre": Scenario("soleil-terre", "la Terre autour du Soleil, un an",
                             soleil_terre, duree=365.25, pas=0.5),
    "solaire": Scenario("solaire", "le Soleil et les huit planètes (Neptune : ~165 ans)",
                        systeme_solaire, duree=4000.0, pas=0.5),
    "huit": Scenario("huit", "trois corps égaux en « huit » (G = 1, période 6,326)",
                     huit, duree=PERIODE_HUIT * 3, pas=0.002, unite=1.0, etiquette="u.t.",
                     unite_longueur=1.0, etiquette_longueur="u.l."),
}
