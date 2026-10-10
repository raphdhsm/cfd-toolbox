"""
Étude de convergence en maillage du cas A (k-ε) de la pièce de Nielsen, pilotée avec foamlib.

Ce que fait le script :
  1. clone le cas cas_kepsilon en trois copies (grossier, moyen, fin) avec foamlib ;
  2. génère le maillage de chaque copie avec maillage_nielsen.py (facteur de raffinement) ;
  3. lance les trois calculs avec foamlib (Allrun : blockMesh, checkMesh, foamRun) ;
  4. extrait les grandeurs d'intérêt avec pyvista et calcule l'indice GCI (Celik et al., 2008).

À placer dans le dossier piece2d, à côté de :
    cas_kepsilon/  maillage_nielsen.py  mesures_nielsen_verticales.csv  mesures_nielsen_horizontales.csv

Utilisation :
    python3 etude_maillage.py            # lance les calculs manquants puis analyse
    python3 etude_maillage.py --analyse  # analyse seulement (calculs déjà faits)

Les cas sont créés dans le dossier etude_maillage/. Un calcul déjà terminé n'est pas relancé.
"""

import os
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from foamlib import FoamCase

# Se placer dans le dossier du script (fonctionne aussi avec le bouton Exécuter de VS Code)
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ".")
from maillage_nielsen import construire_maillage  # noqa: E402  (notre script de maillage)

# =============================================================================
# 1. PARAMÈTRES DE L'ÉTUDE
# =============================================================================
CAS_MODELE = "cas_kepsilon"     # cas de référence à cloner (conditions limites, schémas...)
DOSSIER = Path("etude_maillage")
# Facteur de raffinement de chaque niveau : environ x1,4 en nombre de mailles par direction
# entre deux niveaux (Celik recommande un rapport > 1,3)
NIVEAUX = {"grossier": 0.67, "moyen": 1.0, "fin": 1.5}
ITERATIONS = {"grossier": 3000, "moyen": 3000, "fin": 5000}   # le maillage fin converge plus lentement

H, L = 3.0, 9.0
h = 0.056 * H
U0 = 0.455
Z = 0.05
COULEURS = {"grossier": "#86b6ef", "moyen": "#3987e5", "fin": "#184f95"}
STYLES = {"grossier": ":", "moyen": "--", "fin": "-"}
COULEUR_EXP = "#eb6834"


# =============================================================================
# 2. PRÉPARATION ET LANCEMENT DES CALCULS (foamlib)
# =============================================================================
def preparer_et_lancer(nom, facteur):
    dossier_cas = DOSSIER / f"A_{nom}"
    log = dossier_cas / "log.foamRun"
    if log.exists() and "End" in log.read_text(errors="ignore")[-200:]:
        print(f"[{nom}] déjà calculé, on passe.")
        return

    print(f"[{nom}] clonage de {CAS_MODELE} vers {dossier_cas}")
    if dossier_cas.exists():
        shutil.rmtree(dossier_cas)
    cas = FoamCase(CAS_MODELE).clone(dossier_cas)        # copie propre (sans résultats)
    shutil.rmtree(dossier_cas / "postProcessing", ignore_errors=True)

    # Maillage : même script que d'habitude, avec le facteur de raffinement
    construire_maillage("A", facteur).write(str(dossier_cas / "system" / "blockMeshDict"))

    # Paramètres de calcul modifiés directement depuis Python, sans ouvrir le fichier
    cas.control_dict["endTime"] = ITERATIONS[nom]
    cas.control_dict["writeInterval"] = 500
    cas.control_dict["purgeWrite"] = 0
    (dossier_cas / "cas.foam").touch()

    print(f"[{nom}] calcul en cours ({ITERATIONS[nom]} itérations)...")
    cas.run()                                              # exécute Allrun
    print(f"[{nom}] terminé.")


# =============================================================================
# 3. POST-TRAITEMENT (pyvista)
# =============================================================================
def lire(dossier_cas, instant=None):
    lecteur = pv.POpenFOAMReader(str(dossier_cas / "cas.foam"))
    lecteur.disable_all_patch_arrays()
    lecteur.enable_patch_array("internalMesh")
    temps = [t for t in lecteur.time_values if t > 0]
    lecteur.set_active_time_value(instant if instant is not None else temps[-1])
    return lecteur.read()["internalMesh"], temps


def profil(maillage, debut, fin, axe):
    ligne = maillage.sample_over_line(debut, fin, resolution=399)
    return ligne.points[:, axe] / H, ligne.point_data["U"][:, 0] / U0


LIGNES = {
    "x = H":   ((H, 0, Z), (H, H, Z), 1),
    "x = 2H":  ((2 * H, 0, Z), (2 * H, H, Z), 1),
    "plafond": ((0, H - h / 2, Z), (L, H - h / 2, Z), 0),
    "sol":     ((0, h / 2, Z), (L, h / 2, Z), 0),
}


def lire_mesures():
    mesures = {}
    for f in ("mesures_nielsen_verticales.csv", "mesures_nielsen_horizontales.csv"):
        lignes = [l.strip() for l in open(f, encoding="utf-8") if l.strip() and not l.startswith("#")][1:]
        for l in lignes:
            nom, c, u, _ = l.split(",")
            mesures.setdefault(nom, []).append((float(c), float(u)))
    return {k: np.array(v).T for k, v in mesures.items()}


def grandeurs(maillage, mesures):
    """Les grandeurs d'intérêt suivies d'un maillage à l'autre."""
    res = {"mailles": maillage.n_cells}
    _, u = profil(maillage, *LIGNES["x = H"])
    res["u_max jet (x=H)"] = u.max()
    _, u = profil(maillage, *LIGNES["x = 2H"])
    res["u_min sol (x=2H)"] = u.min()
    for nom, (debut, fin, axe) in LIGNES.items():
        c, u = profil(maillage, debut, fin, axe)
        c_exp, u_exp = mesures[nom]
        res[f"écart moyen {nom}"] = np.mean(np.abs(np.interp(c_exp, c, u) - u_exp))
    return res


def gci(phi1, phi2, phi3, N1, N2, N3):
    """Indice de convergence en maillage (Celik et al. 2008, J. Fluids Eng. 130).

    1 = maillage fin, 2 = moyen, 3 = grossier ; N = nombre de mailles (cas 2D).
    Renvoie (ordre observé p, valeur extrapolée, GCI du maillage fin en %).
    """
    r21, r32 = (N1 / N2) ** 0.5, (N2 / N3) ** 0.5          # 2D : racine carrée
    e21, e32 = phi2 - phi1, phi3 - phi2
    if e21 == 0 or e32 == 0:
        return np.nan, phi1, 0.0
    s = np.sign(e32 / e21)                                   # s < 0 : convergence oscillante
    p = 2.0
    for _ in range(100):                                     # résolution par point fixe
        q = np.log((r21 ** p - s) / (r32 ** p - s))
        p = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
    phi_ext = (r21 ** p * phi1 - phi2) / (r21 ** p - 1)
    gci_fin = 1.25 * abs((phi1 - phi2) / phi1) / (r21 ** p - 1) * 100
    return p, phi_ext, gci_fin


def analyser():
    mesures = lire_mesures()
    resultats = {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.5), sharey=True)

    for nom in NIVEAUX:
        dossier_cas = DOSSIER / f"A_{nom}"
        maillage, temps = lire(dossier_cas)
        resultats[nom] = grandeurs(maillage, mesures)

        # Contrôle de convergence itérative : écart entre les deux dernières écritures
        if len(temps) >= 2:
            avant, _ = lire(dossier_cas, temps[-2])
            _, u1 = profil(maillage, *LIGNES["x = H"])
            _, u0 = profil(avant, *LIGNES["x = H"])
            resultats[nom]["écart itératif"] = np.max(np.abs(u1 - u0))

        for ax, ligne in zip(axes, ("x = H", "x = 2H")):
            y, u = profil(maillage, *LIGNES[ligne])
            ax.plot(u, y, color=COULEURS[nom], linestyle=STYLES[nom], linewidth=2,
                    label=f"{nom} ({resultats[nom]['mailles']} mailles)")

    for ax, ligne in zip(axes, ("x = H", "x = 2H")):
        y_exp, u_exp = mesures[ligne]
        ax.plot(u_exp, y_exp, "o", color=COULEUR_EXP, markersize=5, label="Mesures (Nielsen 1990)")
        ax.set_title(f"Verticale {ligne}")
        ax.set_xlabel("u / U0")
        ax.axvline(0, color="#999999", linewidth=0.8)
        ax.grid(True, color="#e5e5e5", linewidth=0.6)
        ax.set_ylim(0, 1)
        for cote in ("top", "right"):
            ax.spines[cote].set_visible(False)
    axes[0].set_ylabel("y / H")
    axes[0].legend(frameon=False, loc="best")
    fig.suptitle("Pièce de Nielsen, cas A (k-ε) : convergence en maillage", fontsize=13)
    fig.tight_layout()
    fig.savefig("etude_maillage.png", dpi=150)

    # ---- Tableau des résultats
    cles = [k for k in resultats["moyen"] if k != "mailles"]
    print("\n" + f"{'grandeur':<26}" + "".join(f"{n:>12}" for n in NIVEAUX))
    print(f"{'mailles':<26}" + "".join(f"{resultats[n]['mailles']:>12d}" for n in NIVEAUX))
    for k in cles:
        print(f"{k:<26}" + "".join(f"{resultats[n].get(k, np.nan):>12.4f}" for n in NIVEAUX))

    # ---- GCI sur les deux grandeurs physiques
    N = [resultats[n]["mailles"] for n in ("fin", "moyen", "grossier")]
    print("\nIndice de convergence en maillage (GCI, maillage fin) :")
    for k in ("u_max jet (x=H)", "u_min sol (x=2H)"):
        phi = [resultats[n][k] for n in ("fin", "moyen", "grossier")]
        p, ext, g = gci(*phi, *N)
        print(f"  {k:<20} ordre observé p = {p:5.2f}   valeur extrapolée = {ext:7.4f}   GCI = {g:5.2f} %")
    print("\nGraphique écrit : etude_maillage.png")


if __name__ == "__main__":
    if "--analyse" not in sys.argv:
        DOSSIER.mkdir(exist_ok=True)
        for nom, facteur in NIVEAUX.items():
            preparer_et_lancer(nom, facteur)
    analyser()
