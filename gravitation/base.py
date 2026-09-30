"""Base de données des corps célestes (partie IV du sujet), en SQLite.

Deux tables : CORPS(id_corps, nom, masse en kg) et
ETAT(id_corps, datem, x, y, z, vx, vy, vz) avec les positions en unités astronomiques,
les vitesses en km/s et `datem` en jours. Un corps peut avoir plusieurs états datés.
"""

import math
import sqlite3

import numpy as np

from .constantes import JOUR, KM, UA
from .systemes import Systeme

SCHEMA = """
CREATE TABLE IF NOT EXISTS CORPS (
    id_corps INTEGER PRIMARY KEY,
    nom      TEXT NOT NULL,
    masse    REAL NOT NULL CHECK (masse > 0)
);
CREATE TABLE IF NOT EXISTS ETAT (
    id_corps INTEGER NOT NULL REFERENCES CORPS(id_corps),
    datem    REAL NOT NULL,
    x REAL NOT NULL, y REAL NOT NULL, z REAL NOT NULL,
    vx REAL NOT NULL, vy REAL NOT NULL, vz REAL NOT NULL,
    PRIMARY KEY (id_corps, datem)
);
"""

# IV.A — la liste des masses
REQUETE_MASSES = "SELECT masse FROM CORPS ORDER BY id_corps"

# IV.B.1 — nombre de corps ayant au moins un état antérieur à tmin
REQUETE_NB_CONNUS = "SELECT COUNT(DISTINCT id_corps) FROM ETAT WHERE datem < :tmin"

# IV.B.2 — pour chaque corps, la date de son dernier état antérieur à tmin
REQUETE_DERNIERE_DATE = """
SELECT id_corps, MAX(datem) AS date_der
FROM ETAT
WHERE datem < :tmin
GROUP BY id_corps
"""

# IV.B.3 — masse, position et vitesse des corps assez massifs situés dans le cube d'arête
# `arete` centré sur l'origine, triés par distance croissante à l'origine
# (le tri sur le carré de la distance donne le même ordre, sans racine carrée).
REQUETE_SELECTION = f"""
SELECT c.id_corps, c.nom, c.masse, e.x, e.y, e.z, e.vx, e.vy, e.vz
FROM CORPS c
JOIN ({REQUETE_DERNIERE_DATE}) d ON d.id_corps = c.id_corps
JOIN ETAT e ON e.id_corps = d.id_corps AND e.datem = d.date_der
WHERE c.masse >= :masse_min
  AND ABS(e.x) <= :arete / 2 AND ABS(e.y) <= :arete / 2 AND ABS(e.z) <= :arete / 2
ORDER BY e.x * e.x + e.y * e.y + e.z * e.z
"""


class Base:
    """Connexion SQLite (en mémoire par défaut) et requêtes du sujet."""

    def __init__(self, chemin=":memory:"):
        self.conn = sqlite3.connect(chemin)
        self.conn.executescript(SCHEMA)

    def fermer(self):
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.fermer()

    # -- écriture --------------------------------------------------------
    def ajouter_corps(self, nom, masse):
        with self.conn:
            cur = self.conn.execute("INSERT INTO CORPS(nom, masse) VALUES (?, ?)", (nom, masse))
        return cur.lastrowid

    def ajouter_etat(self, id_corps, datem, position_m, vitesse_ms):
        """Enregistre un état donné en SI (m, m/s) ; la base le stocke en UA et km/s."""
        ligne = (int(id_corps), float(datem),
                 *map(float, np.asarray(position_m) / UA), *map(float, np.asarray(vitesse_ms) / KM))
        with self.conn:
            self.conn.execute("INSERT OR REPLACE INTO ETAT VALUES (?,?,?,?,?,?,?,?)", ligne)

    def enregistrer_systeme(self, systeme, datem=0.0):
        """Ajoute chaque corps du système et son état à la date `datem` (jours).
        Retourne les identifiants créés."""
        ids = []
        for nom, m, p, v in zip(systeme.noms, systeme.masses, systeme.positions,
                                systeme.vitesses):
            i = self.ajouter_corps(nom, float(m))
            self.ajouter_etat(i, datem, p, v)
            ids.append(i)
        return ids

    def enregistrer_trajectoire(self, ids, traj, datem0=0.0):
        """Stocke tous les instants d'une trajectoire (t en secondes → datem en jours)."""
        lignes = []
        for t, ps, vs in zip(traj.t, traj.positions, traj.vitesses):
            for i, p, v in zip(ids, ps, vs):
                lignes.append((i, datem0 + float(t) / JOUR, *map(float, p / UA),
                               *map(float, v / KM)))
        with self.conn:
            self.conn.executemany("INSERT OR REPLACE INTO ETAT VALUES (?,?,?,?,?,?,?,?)", lignes)

    # -- requêtes du sujet ----------------------------------------------
    def masses(self):
        return [r[0] for r in self.conn.execute(REQUETE_MASSES)]

    def nb_corps_connus_avant(self, tmin):
        return self.conn.execute(REQUETE_NB_CONNUS, {"tmin": tmin}).fetchone()[0]

    def nb_corps(self):
        return self.conn.execute("SELECT COUNT(*) FROM CORPS").fetchone()[0]

    def dernieres_dates(self, tmin):
        return dict(self.conn.execute(REQUETE_DERNIERE_DATE, {"tmin": tmin}).fetchall())

    def selectionner(self, tmin, masse_min=0.0, arete=math.inf):
        """Le système des corps de masse ≥ masse_min dans le cube d'arête `arete` (UA),
        dans l'état connu le plus récent avant `tmin`, prêt à simuler (unités SI).

        Retourne le système et la liste des identifiants de ses corps, dans le même ordre."""
        arete = 1e300 if math.isinf(arete) else arete
        lignes = self.conn.execute(REQUETE_SELECTION, {"tmin": tmin, "masse_min": masse_min,
                                                       "arete": arete}).fetchall()
        if not lignes:
            raise LookupError("aucun corps ne correspond à la sélection")
        systeme = Systeme([l[1] for l in lignes], [l[2] for l in lignes],
                          [[c * UA for c in l[3:6]] for l in lignes],
                          [[c * KM for c in l[6:9]] for l in lignes])
        return systeme, [l[0] for l in lignes]
