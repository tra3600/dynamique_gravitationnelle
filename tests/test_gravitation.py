import math

import numpy as np
import pytest

from gravitation import (Systeme, accelerations, bilan, energie, etat_suiv, force2, forceN,
                         pos_suiv, simuler, simuler_reference)
from gravitation import bench, diagnostics, edo, systemes
from gravitation.base import Base
from gravitation.cli import main
from gravitation.constantes import G, JOUR, UA
from gravitation.integrateurs import METHODES
from gravitation.vecteurs import produit_scalaire, smul, vdif, vsom


# ---------------------------------------------------------------- vecteurs
def test_vecteurs():
    assert smul(2, [1, 2, 3]) == [2, 4, 6]
    assert vsom([1, 2, 3], [4, 5, 6]) == [5, 7, 9]
    assert vdif([1, 2, 3], [4, 5, 6]) == [-3, -3, -3]
    assert produit_scalaire([1, 2, 3], [4, 5, 6]) == 32
    with pytest.raises(ValueError):
        vsom([1], [1, 2])


# ---------------------------------------------------------------- EDO
def test_euler_energie_croit_du_facteur_du_sujet():
    omega, h_n = 2 * math.pi, 100
    y, z = edo.euler(edo.oscillateur(omega), 3.0, 0.0, 0.0, 3.0, h_n)
    E = edo.energie_oscillateur(y, z, omega)
    h = 3.0 / (h_n - 1)
    for a, b in zip(E, E[1:]):                      # E_{i+1} = (1 + h²ω²)·E_i
        assert b == pytest.approx(a * (1 + (h * omega) ** 2), rel=1e-9)


def test_verlet_conserve_mieux_l_energie_qu_euler():
    omega = 2 * math.pi
    f = edo.oscillateur(omega)
    E_e = edo.energie_oscillateur(*edo.euler(f, 3.0, 0.0, 0.0, 3.0, 100), omega)
    E_v = edo.energie_oscillateur(*edo.verlet(f, 3.0, 0.0, 0.0, 3.0, 100), omega)
    assert max(abs(e / E_v[0] - 1) for e in E_v) < 0.05
    assert E_e[-1] / E_e[0] > 5


def test_edo_valide_les_arguments():
    with pytest.raises(ValueError):
        edo.euler(lambda y: y, 1, 0, 0, 1, 1)
    with pytest.raises(ValueError):
        edo.verlet(lambda y: y, 1, 0, 1, 0, 10)


def test_verlet_accepte_les_tableaux():
    y, z = edo.verlet(lambda y: -y, np.array([1.0, 0.0]), np.array([0.0, 1.0]), 0, 1, 50)
    assert y[-1].shape == (2,)


# ---------------------------------------------------------------- forces
def test_force2_terre_lune():
    F = force2(5.972e24, [0, 0, 0], 7.348e22, [3.844e8, 0, 0])
    assert F[0] == pytest.approx(G * 5.972e24 * 7.348e22 / 3.844e8 ** 2, rel=1e-12)
    assert F[1] == F[2] == 0.0


def test_force2_positions_confondues():
    with pytest.raises(ValueError):
        force2(1, [0, 0, 0], 1, [0, 0, 0])


def test_troisieme_loi_de_newton():
    m, p = [3.0, 5.0, 7.0], [[0, 0, 0], [1, 2, 0], [-2, 1, 3]]
    somme = np.sum([forceN(j, m, p) for j in range(3)], axis=0)
    assert np.allclose(somme, 0, atol=1e-18)


def test_accelerations_coincident_avec_forceN():
    s = systemes.amas(6, graine=3)
    m, p, _ = s.en_listes()
    a = accelerations(s.masses, s.positions)
    for j in range(6):
        assert a[j] == pytest.approx(np.array(forceN(j, m, p)) / m[j], rel=1e-10)


def test_accelerations_adoucissement_borne_la_force():
    m, p = [1.0, 1.0], [[0, 0, 0], [1e-9, 0, 0]]
    assert np.all(np.isfinite(accelerations(m, p, adoucissement=1.0)))
    with pytest.raises(ValueError):
        accelerations(m, [[0, 0, 0], [0, 0, 0]])


# ---------------------------------------------------------------- intégrateurs
def test_etat_suiv_listes_egal_a_un_pas_vectorise():
    s = systemes.terre_lune()
    m, p, v = s.en_listes()
    pos, vit = etat_suiv(m, p, v, 60.0)
    assert np.allclose(pos_suiv(m, p, v, 60.0), pos)
    traj = simuler(s, 60.0, 60.0, "verlet")
    assert np.allclose(traj.positions[-1], pos, rtol=1e-12)
    assert np.allclose(traj.vitesses[-1], vit, rtol=1e-9, atol=1e-9)


def test_methode_inconnue_ou_pas_invalide():
    s = systemes.terre_lune()
    with pytest.raises(ValueError):
        simuler(s, 10, 1, "inconnue")
    with pytest.raises(ValueError):
        simuler(s, 10, 0)


def test_lune_periode_et_conservation():
    s = systemes.terre_lune()
    traj = simuler(s, 100 * JOUR, 0.02 * JOUR, "verlet", enregistrer_tous=5)
    b = bilan(traj)
    assert b["energie"] < 1e-8 and b["moment_cinetique"] < 1e-8 and b["quantite_de_mouvement"] < 1e-8
    a = 3.844e8
    T_theorique = 2 * math.pi * math.sqrt(a ** 3 / (G * (5.972e24 + 7.348e22)))
    rel = traj.positions[:, 1] - traj.positions[:, 0]
    T = diagnostics.periode_moyenne(traj.t, rel[:, 0])
    assert T == pytest.approx(T_theorique, rel=1e-3)


def test_euler_derive_verlet_non():
    s = systemes.terre_lune()
    e = bilan(simuler(s, 30 * JOUR, 0.05 * JOUR, "euler"))["energie"]
    v = bilan(simuler(s, 30 * JOUR, 0.05 * JOUR, "verlet"))["energie"]
    assert e > 1e3 * v


def test_ordre_de_convergence_verlet_et_rk4():
    s = systemes.terre_lune()
    ref = simuler_reference(s, 5 * JOUR, nb_points=2)
    erreur = {}
    for methode in ("verlet", "rk4"):
        for h in (0.1, 0.05):
            t = simuler(s, 5 * JOUR, h * JOUR, methode)
            erreur[methode, h] = np.max(np.abs(t.positions[-1] - ref.positions[-1]))
    assert erreur["verlet", 0.1] / erreur["verlet", 0.05] == pytest.approx(4, rel=0.3)   # ordre 2
    assert erreur["rk4", 0.1] / erreur["rk4", 0.05] > 10                                  # ordre 4


def test_huit_est_periodique():
    s = systemes.huit()
    traj = simuler(s, systemes.PERIODE_HUIT, 1e-3, "rk4")
    # conditions initiales connues à 8 chiffres : retour à ~1e-4 près après une période
    assert np.allclose(traj.positions[-1], s.positions, atol=1e-3)
    assert bilan(traj)["energie"] < 1e-9


def test_solaire_conserve_energie_et_barycentre():
    s = systemes.systeme_solaire()
    assert s.n == 9
    traj = simuler(s, 2000 * JOUR, 0.5 * JOUR, "verlet", enregistrer_tous=20)
    b = bilan(traj)
    assert b["energie"] < 1e-6 and b["moment_cinetique"] < 1e-6
    assert np.allclose(np.linalg.norm(s.barycentre()), 0, atol=1)


def test_enregistrer_tous_garde_le_dernier_instant():
    s = systemes.terre_lune()
    traj = simuler(s, 10 * 60.0, 60.0, "verlet", enregistrer_tous=4)
    assert list(traj.t) == [0, 240, 480, 600]


def test_toutes_les_methodes_avancent():
    s = systemes.terre_lune()
    for methode in METHODES:
        assert not np.allclose(simuler(s, 600, 60, methode).positions[-1], s.positions)


def test_systeme_valide_ses_entrees():
    with pytest.raises(ValueError):
        Systeme(["a"], [1.0], [[0, 0]], [[0, 0, 0]])
    with pytest.raises(ValueError):
        Systeme(["a"], [-1.0], [[0, 0, 0]], [[0, 0, 0]])


def test_referentiel_barycentrique():
    s = systemes.soleil_terre()
    assert np.allclose(s.barycentre(), 0, atol=1e-3)
    assert np.allclose(s.vitesse_barycentre(), 0, atol=1e-9)
    assert energie(s) < 0                                  # système lié


# ---------------------------------------------------------------- base SQL
@pytest.fixture
def base():
    with Base() as b:
        b.enregistrer_systeme(systemes.systeme_solaire(), datem=0.0)
        yield b


def test_requetes_du_sujet(base):
    assert len(base.masses()) == 9
    assert base.nb_corps_connus_avant(1.0) == 9
    assert base.nb_corps_connus_avant(0.0) == 0            # « antérieur » au sens strict
    assert set(base.dernieres_dates(1.0).values()) == {0.0}


def test_selection_filtre_et_ordonne(base):
    s, ids = base.selectionner(1.0, masse_min=1e26, arete=100)
    assert set(s.noms) == {"Soleil", "Jupiter", "Saturne", "Neptune"}      # Uranus : 8,7e25 kg < 1e26
    assert len(ids) == s.n
    d = np.linalg.norm(s.positions, axis=1)
    assert np.all(np.diff(d) >= 0)                          # triés par distance à l'origine
    petit, _ = base.selectionner(1.0, arete=12)             # Uranus et Neptune sont hors du cube
    assert "Neptune" not in petit.noms and "Terre" in petit.noms


def test_selection_prend_le_dernier_etat_connu(base):
    s0, ids0 = base.selectionner(1.0)
    traj = simuler(s0, 100 * JOUR, JOUR, "verlet")
    base.enregistrer_trajectoire(ids0, traj, datem0=0.0)
    s, ids = base.selectionner(50.5)                        # dernier état avant 50,5 j : t = 50
    i = list(traj.t / JOUR).index(50.0)
    par_id = dict(zip(ids, s.positions))
    for k, ident in enumerate(ids0):
        assert np.allclose(par_id[ident], traj.positions[i, k], rtol=1e-9, atol=1e-3)


def test_aller_retour_unites(base):
    s, _ = base.selectionner(1.0)
    ref = systemes.systeme_solaire()
    assert sorted(s.noms) == sorted(ref.noms)
    terre = s.noms.index("Terre"), ref.noms.index("Terre")
    assert np.allclose(s.positions[terre[0]], ref.positions[terre[1]], rtol=1e-12, atol=1)
    d = np.linalg.norm(s.positions[terre[0]] - s.positions[s.noms.index("Soleil")])
    assert d == pytest.approx(UA, rel=1e-9)


def test_selection_vide(base):
    with pytest.raises(LookupError):
        base.selectionner(1.0, masse_min=1e40)


# ---------------------------------------------------------------- complexité / CLI
def test_complexite_quadratique_pour_la_version_listes():
    tailles = [10, 20, 40, 80]
    t = bench.mesurer(tailles, repetitions=2)
    assert 1.5 < bench.pente_loglog(tailles, t["listes"]) < 2.5


@pytest.mark.parametrize("argv", [
    ["oscillateur"],
    ["orbite", "terre-lune", "--duree", "5"],
    ["comparer", "terre-lune", "--duree", "5"],
    ["comparer", "huit", "--sans-reference", "--duree", "3"],
    ["cas", "--methode", "rk4", "--duree", "3"],
    ["base", "--jours", "30"],
])
def test_cli(argv, capsys):
    assert main(argv) == 0
    assert capsys.readouterr().out


def test_cli_graphiques(tmp_path):
    pytest.importorskip("matplotlib")
    assert main(["orbite", "terre-lune", "--duree", "5", "--sauver", str(tmp_path)]) == 0
    assert main(["oscillateur", "--sauver", str(tmp_path)]) == 0
    assert (tmp_path / "orbite_terre-lune.png").exists() and (tmp_path / "oscillateur.png").exists()
