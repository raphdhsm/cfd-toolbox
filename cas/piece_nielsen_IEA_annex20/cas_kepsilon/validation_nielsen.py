"""
Validation de la pièce de Nielsen : simulation OpenFOAM contre les mesures laser-Doppler.

Compare u/U0 sur les quatre lignes de mesure de l'IEA Annex 20 :
    - deux verticales  : x = H et x = 2H
    - deux horizontales: y = H - h/2 (dans le jet, sous le plafond) et y = h/2 (près du sol)

Utilisation (dans le dossier du cas, là où se trouve cas.foam) :
    python3 validation_nielsen.py            # dernier instant écrit
    python3 validation_nielsen.py 3000       # un instant précis

Fichiers nécessaires dans le même dossier :
    mesures_nielsen_verticales.csv, mesures_nielsen_horizontales.csv

Sorties :
    validation_nielsen.png   le graphique (4 panneaux)
    écarts simulation/mesure affichés dans le terminal
"""

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
NOM_CAS = "cas A (k-ε, lois de paroi)"

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
    "plafond": "Horizontale y = H − h/2 (dans le jet)",
    "sol": "Horizontale y = h/2 (près du sol)",
}

COULEUR_SIM = "#2a78d6"   # simulation : trait continu bleu
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


def main():
    mesures = lire_mesures()

    lecteur = pv.POpenFOAMReader("cas.foam")
    lecteur.disable_all_patch_arrays()
    lecteur.enable_patch_array("internalMesh")
    instant = float(sys.argv[1]) if len(sys.argv) > 1 else max(lecteur.time_values)
    lecteur.set_active_time_value(instant)
    maillage = lecteur.read()["internalMesh"]
    print(f"Instant analysé : {instant:g}")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    print("\nÉcart simulation - mesure sur u/U0 (aux points de mesure) :")
    print(f"  {'ligne':<10} {'écart moyen':>12} {'écart max':>10}")

    for ax, (nom, debut, fin, sens) in zip(axes.flat, LIGNES):
        coord_sim, u_sim = echantillonner(maillage, debut, fin, sens)
        coord_exp, u_exp, _ = mesures[nom]

        # Écart : on interpole la simulation aux positions des mesures
        u_sim_aux_mesures = np.interp(coord_exp, coord_sim, u_sim)
        ecart = np.abs(u_sim_aux_mesures - u_exp)
        print(f"  {nom:<10} {ecart.mean():>12.3f} {ecart.max():>10.3f}")

        if sens == "vertical":
            ax.plot(u_sim, coord_sim, color=COULEUR_SIM, linewidth=2, label="OpenFOAM")
            ax.plot(u_exp, coord_exp, "o", color=COULEUR_EXP, markersize=6,
                    label="Mesures (Nielsen 1990)")
            ax.set_xlabel("u / U0")
            ax.set_ylabel("y / H")
            ax.set_ylim(0, 1)
            ax.axvline(0, color="#999999", linewidth=0.8)
        else:
            ax.plot(coord_sim, u_sim, color=COULEUR_SIM, linewidth=2, label="OpenFOAM")
            ax.plot(coord_exp, u_exp, "o", color=COULEUR_EXP, markersize=6,
                    label="Mesures (Nielsen 1990)")
            ax.set_xlabel("x / H")
            ax.set_ylabel("u / U0")
            ax.set_xlim(0, 3)
            ax.axhline(0, color="#999999", linewidth=0.8)

        ax.set_title(TITRES[nom], fontsize=12)
        ax.grid(True, color="#e5e5e5", linewidth=0.6)
        for cote in ("top", "right"):
            ax.spines[cote].set_visible(False)

    axes.flat[0].legend(frameon=False, loc="best")
    fig.suptitle(f"Pièce de Nielsen (IEA Annex 20), {NOM_CAS} : simulation contre mesures",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig("validation_nielsen.png", dpi=150)
    print("\nGraphique écrit : validation_nielsen.png")


if __name__ == "__main__":
    main()
