"""Interface en ligne de commande : `python -m gravitation <commande>`.

    python -m gravitation cas                               # tous les scénarios, bilan des invariants
    python -m gravitation orbite solaire --graphique
    python -m gravitation comparer huit --sauver figures    # Euler, Verlet, RK4 vs référence scipy
    python -m gravitation oscillateur --graphique           # II.B.3 : Euler diverge, Verlet non
    python -m gravitation complexite --nmax 200             # III.B.4-5 : pente log-log
    python -m gravitation base                              # partie IV : base SQL de corps
"""

import argparse
import os
import sys

import numpy as np

from . import __version__, bench, diagnostics, edo
from .base import Base
from .constantes import JOUR
from .integrateurs import METHODES, simuler, simuler_reference
from .systemes import SCENARIOS


def _montrer(fig, sauver, nom):
    import matplotlib.pyplot as plt
    if sauver:
        os.makedirs(sauver, exist_ok=True)
        chemin = os.path.join(sauver, f"{nom}.png")
        fig.savefig(chemin, dpi=110)
        print(f"   🖼  {chemin}")
        plt.close(fig)
    else:
        plt.show()


def _preparer_matplotlib(args):
    if getattr(args, "sauver", None):
        import matplotlib
        matplotlib.use("Agg")


def _reglages(args, scenario):
    """(durée, pas) en secondes ; --duree et --pas sont exprimés dans l'unité du scénario."""
    duree = (args.duree if args.duree is not None else scenario.duree) * scenario.unite
    h = (args.pas if args.pas is not None else scenario.pas) * scenario.unite
    n = max(1, round(duree / h))
    return n * h, h


def _garder(traj_longueur, max_points=2000):
    return max(1, traj_longueur // max_points)


def _simuler(scenario, systeme, duree, h, methode):
    n = round(duree / h)
    return simuler(systeme, duree, h, methode, enregistrer_tous=_garder(n))


def _afficher_bilan(traj):
    b = diagnostics.bilan(traj)
    print(f"   {traj.methode:<20} dérive E {b['energie']:>9.2e}   L {b['moment_cinetique']:>9.2e}"
          f"   p {b['quantite_de_mouvement']:>9.2e}   {traj.duree_calcul:>7.2f} s")


def _figure_orbites(traj, scenario, titre):
    from . import visualisation
    return visualisation.tracer_orbites(traj, scenario.unite_longueur,
                                        scenario.etiquette_longueur, titre)


# ----------------------------------------------------------------------
def cmd_cas(args):
    """Tous les scénarios avec le schéma choisi : quelle précision, en combien de temps ?"""
    _preparer_matplotlib(args)
    print("═" * 78)
    print(f"  CAS ILLUSTRÉS — schéma {args.methode}")
    print("═" * 78)
    for nom, sc in SCENARIOS.items():
        duree, h = _reglages(args, sc)
        systeme = sc.construire()
        print(f"\n── {nom} : {sc.description}")
        print(f"   {systeme.n} corps, {round(duree / h)} pas de {h / sc.unite:g} {sc.etiquette}")
        traj = _simuler(sc, systeme, duree, h, args.methode)
        _afficher_bilan(traj)
        if args.graphique or args.sauver:
            _montrer(_figure_orbites(traj, sc, f"{nom} — {sc.description}"), args.sauver,
                     f"orbite_{nom}")


def cmd_orbite(args):
    _preparer_matplotlib(args)
    sc = SCENARIOS[args.scenario]
    duree, h = _reglages(args, sc)
    traj = _simuler(sc, sc.construire(), duree, h, args.methode)
    print(f"{sc.description} — {duree / sc.unite:g} {sc.etiquette}, "
          f"pas {h / sc.unite:g} {sc.etiquette}")
    _afficher_bilan(traj)
    if args.graphique or args.sauver:
        _montrer(_figure_orbites(traj, sc, sc.description), args.sauver, f"orbite_{args.scenario}")


def cmd_comparer(args):
    """Euler, Euler symplectique, Verlet et RK4 face à la solution de référence scipy."""
    _preparer_matplotlib(args)
    from . import visualisation
    sc = SCENARIOS[args.scenario]
    duree, h = _reglages(args, sc)
    systeme = sc.construire()
    echelle = float(np.max(np.linalg.norm(systeme.positions, axis=1)))
    print(f"{sc.description} — {duree / sc.unite:g} {sc.etiquette}, "
          f"pas {h / sc.unite:g} {sc.etiquette}")
    ref = None
    if not args.sans_reference:
        ref = simuler_reference(systeme, duree, nb_points=2)
        print(f"   {'référence DOP853':<20} {'':>36} {ref.duree_calcul:>7.2f} s")
    trajs = []
    for methode in METHODES:
        traj = _simuler(sc, systeme, duree, h, methode)
        trajs.append(traj)
        b = diagnostics.bilan(traj)
        err = ""
        if ref is not None:
            e = np.max(np.linalg.norm(traj.positions[-1] - ref.positions[-1], axis=1)) / echelle
            err = f"   erreur position finale {e:>9.2e}"
        print(f"   {methode:<20} dérive E {b['energie']:>9.2e}{err}   {traj.duree_calcul:>7.2f} s")
    if args.graphique or args.sauver:
        _montrer(visualisation.tracer_energie(trajs, sc.unite, sc.etiquette,
                                              f"Énergie — {args.scenario}"),
                 args.sauver, f"energie_{args.scenario}")


def cmd_oscillateur(args):
    _preparer_matplotlib(args)
    omega = args.omega
    f = edo.oscillateur(omega)
    y_e, z_e = edo.euler(f, args.y0, 0.0, 0.0, args.tmax, args.n)
    y_v, z_v = edo.verlet(f, args.y0, 0.0, 0.0, args.tmax, args.n)
    h = args.tmax / (args.n - 1)
    print(f"y'' = −ω²y, ω = {omega:.4f}, y0 = {args.y0}, h = {h:.4f}, n = {args.n}")
    for nom, y, z in (("Euler", y_e, z_e), ("Verlet", y_v, z_v)):
        E = edo.energie_oscillateur(y, z, omega)
        print(f"   {nom:<8} E_final / E_0 = {E[-1] / E[0]:.4f}   "
              f"(max |E/E_0 − 1| = {max(abs(e / E[0] - 1) for e in E):.2e})")
    print(f"   Euler : E_(i+1) = (1 + h²ω²)·E_i, soit un facteur {(1 + (h * omega) ** 2) ** (args.n - 1):.4f} "
          f"sur {args.n - 1} pas.")
    if args.graphique or args.sauver:
        from . import visualisation
        t = edo.grille(0.0, args.tmax, args.n)
        exact = [args.y0 * np.cos(omega * ti) for ti in t]
        _montrer(visualisation.tracer_oscillateur(t, (y_e, z_e), (y_v, z_v), exact, omega),
                 args.sauver, "oscillateur")


def cmd_complexite(args):
    _preparer_matplotlib(args)
    tailles = np.unique(np.geomspace(args.nmin, args.nmax, args.points).astype(int))
    mesures = bench.mesurer(tailles, repetitions=args.repetitions)
    pentes = {k: bench.pente_loglog(tailles, v) for k, v in mesures.items()}
    print(f"{'N':>6}{'listes (sujet)':>18}{'numpy':>14}{'gain':>9}")
    for i, n in enumerate(tailles):
        a, b = mesures["listes"][i], mesures["numpy"][i]
        print(f"{n:>6}{a:>17.2e}s{b:>13.2e}s{a / b:>8.0f}×")
    print(f"\nPentes log-log : listes {pentes['listes']:.2f}, numpy {pentes['numpy']:.2f} "
          "(N·(N−1) interactions par pas : O(N²) ; "
          "la pente numpy est plus faible tant que le surcoût fixe domine, aux petits N)")
    if args.graphique or args.sauver:
        from . import visualisation
        _montrer(visualisation.tracer_complexite(tailles, mesures, pentes), args.sauver,
                 "complexite")


def cmd_base(args):
    """Partie IV : remplit une base SQL, l'interroge, simule la sélection et réécrit les états."""
    _preparer_matplotlib(args)
    sc = SCENARIOS["solaire"]
    with Base(args.fichier or ":memory:") as base:
        if base.nb_corps() == 0:
            base.enregistrer_systeme(sc.construire(), datem=0.0)
        tmin = args.tmin
        print(f"{base.nb_corps()} corps ; masses (kg) : "
              f"{', '.join(f'{m:.2e}' for m in base.masses())}")
        print(f"Corps avec un état antérieur à t = {tmin} j : {base.nb_corps_connus_avant(tmin)}")
        print(f"Dernières dates connues : {base.dernieres_dates(tmin)}")
        systeme, ids = base.selectionner(tmin, masse_min=args.masse_min, arete=args.arete)
        print(f"Sélection (masse ≥ {args.masse_min:g} kg, cube de {args.arete:g} UA) : "
              f"{', '.join(systeme.noms)}")
        duree = args.jours * JOUR
        traj = simuler(systeme, duree, 0.5 * JOUR, "verlet", enregistrer_tous=10)
        _afficher_bilan(traj)
        base.enregistrer_trajectoire(ids, traj, datem0=tmin)
        n_etats = base.conn.execute("SELECT COUNT(*) FROM ETAT").fetchone()[0]
        print(f"Trajectoire enregistrée : {n_etats} états dans la base.")
        if args.graphique or args.sauver:
            _montrer(_figure_orbites(traj, sc, "Simulation depuis la base"), args.sauver,
                     "base_orbites")


# ----------------------------------------------------------------------
def _parser():
    p = argparse.ArgumentParser(prog="gravitation", description=__doc__.split("\n")[0])
    p.add_argument("--version", action="version", version=f"gravitation {__version__}")
    sub = p.add_subparsers(dest="commande", required=True)

    def commun(sp, methode=True):
        sp.add_argument("--graphique", action="store_true", help="affiche les graphiques")
        sp.add_argument("--sauver", metavar="DOSSIER", help="enregistre les graphiques en PNG")
        if methode:
            sp.add_argument("--methode", choices=list(METHODES), default="verlet")
        return sp

    def reglages(sp):
        sp.add_argument("--duree", type=float, help="durée (dans l'unité du scénario : jours…)")
        sp.add_argument("--pas", type=float, help="pas de temps (même unité)")

    sp = commun(sub.add_parser("cas", help="tous les scénarios, bilan des invariants"))
    reglages(sp)
    sp.set_defaults(fonction=cmd_cas)

    sp = commun(sub.add_parser("orbite", help="simule un scénario"))
    sp.add_argument("scenario", choices=list(SCENARIOS))
    reglages(sp)
    sp.set_defaults(fonction=cmd_orbite)

    sp = commun(sub.add_parser("comparer", help="schémas d'intégration vs référence scipy"),
                methode=False)
    sp.add_argument("scenario", choices=list(SCENARIOS))
    sp.add_argument("--sans-reference", action="store_true", help="n'appelle pas scipy")
    reglages(sp)
    sp.set_defaults(fonction=cmd_comparer)

    sp = commun(sub.add_parser("oscillateur", help="Euler contre Verlet sur y'' = −ω²y"),
                methode=False)
    sp.add_argument("--omega", type=float, default=2 * np.pi)
    sp.add_argument("--y0", type=float, default=3.0)
    sp.add_argument("--tmax", type=float, default=3.0)
    sp.add_argument("--n", type=int, default=100)
    sp.set_defaults(fonction=cmd_oscillateur)

    sp = commun(sub.add_parser("complexite", help="temps d'un pas en fonction de N"),
                methode=False)
    sp.add_argument("--nmin", type=int, default=5)
    sp.add_argument("--nmax", type=int, default=200)
    sp.add_argument("--points", type=int, default=8)
    sp.add_argument("--repetitions", type=int, default=3)
    sp.set_defaults(fonction=cmd_complexite)

    sp = commun(sub.add_parser("base", help="base SQL de corps célestes"), methode=False)
    sp.add_argument("--fichier", help="fichier SQLite (en mémoire par défaut)")
    sp.add_argument("--tmin", type=float, default=1.0, help="date de départ (jours)")
    sp.add_argument("--masse-min", type=float, default=1e23, help="masse minimale (kg)")
    sp.add_argument("--arete", type=float, default=80.0, help="arête du cube (UA)")
    sp.add_argument("--jours", type=float, default=365.25)
    sp.set_defaults(fonction=cmd_base)
    return p


def main(argv=None):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = _parser().parse_args(argv)
    try:
        args.fonction(args)
    except (ValueError, LookupError) as e:
        print(f"Erreur : {e}", file=sys.stderr)
        return 1
    return 0
