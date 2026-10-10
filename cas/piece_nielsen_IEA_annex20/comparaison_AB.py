"""
Comparaison des cas A (k-epsilon) et B (k-omega SST) avec les mesures de Nielsen.

A placer dans le dossier piece2d (celui qui contient cas_kepsilon et cas_komegaSST),
avec les deux fichiers de mesures CSV.

Utilisation :
    python3 comparaison_AB.py

Sorties :
    comparaison_AB.png   4 panneaux : mesures + cas A + cas B
    tableau des écarts simulation/mesure pour chaque cas dans le terminal
"""

import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# =============================================================================
# 1. PARAMÈTRES DU CAS
# =============================================================================
H = 3.0
L = 3.0 * H
h = 0.056 * H
U0 = 0.455
Z_MILIEU = 0.05
N_POINTS = 400
# Cas à comparer : (dossier, légende, couleur, style de trait)
CAS = [
    ("cas_kepsilon",  "A : k-ε + lois de paroi", "#2a78d6", "--"),
    ("cas_komegaSST", "B : k-ω SST, y+ ≈ 1",    "#1baf7a", "-"),
]

# Lignes de mesure : (nom, point de départ, point d'arrivée, axe tracé)
LIGNES = [
    ("x = H",       (H, 0, Z_MILIEU),      (H, H, Z_MILIEU),      "vertical"),
    ("x = 2H",      (2 * H, 0, Z_MILIEU),  (2 * H, H, Z_MILIEU),  "vertical"),
    ("plafond",     (0, H - h / 2, Z_MILIEU), (L, H - h / 2, Z_MILIEU), "horizontal"),
    ("sol",         (0, h / 2, Z_MILIEU),  (L, h / 2, Z_MILIEU),  "horizontal"),
]
TITRES = {
    "x = H": "Verticale x = H",
    "x = 2H": "Verticale x = 2H",
    "plafond": "Horizontale y = H - h/2 (dans le jet)",
    "sol": "Horizontale y = h/2 (près du sol)",
}

COULEUR_EXP = "#eb6834"   # mesures : points orange


def lire_mesures():
    """Charge les deux fichiers CSV de mesures dans un dictionnaire {ligne: (coord, u, rms)}."""
    mesures = {}
    for fichier in ("mesures_nielsen_verticales.csv", "mesures_nielsen_horizontales.csv"):
        if not Path(fichier).exists():
            sys.exit(f"Fichier manquant : {fichier} (à mettre dans le dossier du cas)")
        # On saute les lignes de commentaire (#) et l'en-tête, puis on range par ligne de mesure
        lignes = [l.strip() for l in open(fichier, encoding="utf-8")
                  if l.strip() and not l.startswith("#")][1:]
        for l in lignes:
            nom, coord, u, rms = l.split(",")
            mesures.setdefault(nom, []).append((float(coord), float(u), float(rms)))
    # Conversion en tableaux numpy : (coordonnée/H, u/U0, rms/U0)
    return {nom: tuple(np.array(v).T) for nom, v in mesures.items()}


def echantillonner(maillage, debut, fin, sens):
    """Vitesse u/U0 le long d'une ligne ; renvoie (coordonnée/H, u/U0)."""
    ligne = maillage.sample_over_line(debut, fin, resolution=N_POINTS - 1)
    axe = 1 if sens == "vertical" else 0
    return ligne.points[:, axe] / H, ligne.point_data["U"][:, 0] / U0


def lire_dernier_instant(dossier):
    """Charge le dernier instant écrit d'un cas OpenFOAM."""
    lecteur = pv.POpenFOAMReader(os.path.join(dossier, "cas.foam"))
    lecteur.disable_all_patch_arrays()
    lecteur.enable_patch_array("internalMesh")
    instant = max(lecteur.time_values)
    lecteur.set_active_time_value(instant)
    return lecteur.read()["internalMesh"], instant


def main():
    # Se placer dans le dossier du script (fonctionne aussi avec le bouton Exécuter de VS Code)
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    mesures = lire_mesures()

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    resultats = {}

    for dossier, legende, couleur, style in CAS:
        if not os.path.exists(os.path.join(dossier, "cas.foam")):
            sys.exit(f"Introuvable : {dossier}/cas.foam")
        maillage, instant = lire_dernier_instant(dossier)
        print(f"{legende} : instant {instant:g}")
        resultats[legende] = {}

        for ax, (nom, debut, fin, sens) in zip(axes.flat, LIGNES):
            coord_sim, u_sim = echantillonner(maillage, debut, fin, sens)
            coord_exp, u_exp, _ = mesures[nom]
            ecart = np.abs(np.interp(coord_exp, coord_sim, u_sim) - u_exp)
            resultats[legende][nom] = (ecart.mean(), ecart.max())
            if sens == "vertical":
                ax.plot(u_sim, coord_sim, color=couleur, linestyle=style, linewidth=2, label=legende)
            else:
                ax.plot(coord_sim, u_sim, color=couleur, linestyle=style, linewidth=2, label=legende)

    # Mesures et mise en forme
    for ax, (nom, debut, fin, sens) in zip(axes.flat, LIGNES):
        coord_exp, u_exp, _ = mesures[nom]
        if sens == "vertical":
            ax.plot(u_exp, coord_exp, "o", color=COULEUR_EXP, markersize=6, label="Mesures (Nielsen 1990)")
            ax.set_xlabel("u / U0"); ax.set_ylabel("y / H"); ax.set_ylim(0, 1)
            ax.axvline(0, color="#999999", linewidth=0.8)
        else:
            ax.plot(coord_exp, u_exp, "o", color=COULEUR_EXP, markersize=6, label="Mesures (Nielsen 1990)")
            ax.set_xlabel("x / H"); ax.set_ylabel("u / U0"); ax.set_xlim(0, 3)
            ax.axhline(0, color="#999999", linewidth=0.8)
        ax.set_title(TITRES[nom], fontsize=12)
        ax.grid(True, color="#e5e5e5", linewidth=0.6)
        for cote in ("top", "right"):
            ax.spines[cote].set_visible(False)

    axes.flat[0].legend(frameon=False, loc="best")
    fig.suptitle("Pièce de Nielsen (IEA Annex 20) : k-ε contre k-ω SST contre mesures", fontsize=13)
    fig.tight_layout()
    fig.savefig("comparaison_AB.png", dpi=150)

    # Tableau des écarts
    print("\nÉcart moyen | max de u/U0 par rapport aux mesures :")
    noms = [l[0] for l in LIGNES]
    print(f"  {'cas':<26}" + "".join(f"{n:>17}" for n in noms))
    for legende, res in resultats.items():
        print(f"  {legende:<26}" + "".join(f"{res[n][0]:>8.3f} |{res[n][1]:>6.3f} " for n in noms))
    print("\nGraphique écrit : comparaison_AB.png")


if __name__ == "__main__":
    main()
