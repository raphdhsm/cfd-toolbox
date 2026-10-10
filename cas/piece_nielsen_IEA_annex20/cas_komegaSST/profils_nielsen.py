"""
Post-traitement de la pièce de Nielsen avec pyvista.

Extrait les profils de vitesse horizontale u/U0 en fonction de y/H sur les deux
verticales de mesure de Nielsen (x = H et x = 2H), pour plusieurs instants du calcul,
et les superpose sur un graphique. Si les courbes se confondent, la solution ne bouge
plus : la convergence est prouvée.

Utilisation (dans le dossier du cas, là où se trouve cas.foam) :
    python3 profils_nielsen.py                  # les 3 derniers instants écrits
    python3 profils_nielsen.py 2000 2500 3000   # des instants précis

Sorties :
    profils_convergence.png   le graphique
    profils_<instant>.csv     les profils bruts (pour Excel / Google Sheets)
"""

import sys

import matplotlib
matplotlib.use("Agg")  # pas de fenêtre : on écrit directement le fichier PNG
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

# =============================================================================
# 1. PARAMÈTRES DU CAS
# =============================================================================
H = 3.0          # hauteur de la pièce (m)
U0 = 0.455       # vitesse d'entrée (m/s), pour adimensionner
Z_MILIEU = 0.05  # le cas 2D fait 0,1 m d'épaisseur : on coupe au milieu
N_POINTS = 300   # nombre de points le long de chaque ligne de mesure

# Lignes de mesure de Nielsen : verticales à x = H et x = 2H
LIGNES = {"x = H": H, "x = 2H": 2 * H}

# Couleurs : une seule teinte, du clair (instant ancien) au foncé (instant récent),
# + des styles de trait différents pour rester lisible en noir et blanc.
COULEURS = ["#86b6ef", "#3987e5", "#184f95", "#0d366b"]
STYLES = [":", "--", "-", "-."]


def lire_cas(fichier="cas.foam"):
    """Ouvre le cas OpenFOAM avec le lecteur de pyvista."""
    lecteur = pv.POpenFOAMReader(fichier)
    lecteur.disable_all_patch_arrays()          # on ne lit que le volume (plus rapide)
    lecteur.enable_patch_array("internalMesh")
    return lecteur


def profil_vertical(maillage, x):
    """Échantillonne U le long de la verticale x, du sol au plafond.

    Renvoie (y/H, u/U0) : y adimensionné et composante horizontale de la vitesse.
    """
    debut = (x, 0.0, Z_MILIEU)
    fin = (x, H, Z_MILIEU)
    ligne = maillage.sample_over_line(debut, fin, resolution=N_POINTS - 1)
    y = ligne.points[:, 1]
    u = ligne.point_data["U"][:, 0]   # composante x de la vitesse
    return y / H, u / U0


def main():
    lecteur = lire_cas()

    # Instants disponibles (on ignore t = 0, l'état initial)
    disponibles = [t for t in lecteur.time_values if t > 0]
    if len(sys.argv) > 1:
        instants = [float(a) for a in sys.argv[1:]]
        manquants = [t for t in instants if t not in disponibles]
        if manquants:
            sys.exit(f"Instants introuvables : {manquants}. Disponibles : {disponibles}")
    else:
        instants = disponibles[-3:]
    instants = instants[-len(COULEURS):]  # au plus 4 courbes par graphique, pour rester lisible
    print(f"Instants analysés : {instants}")

    # Un graphique par ligne de mesure, côte à côte, même axe y/H
    fig, axes = plt.subplots(1, len(LIGNES), figsize=(10, 5.5), sharey=True)

    for i, t in enumerate(instants):
        lecteur.set_active_time_value(t)
        maillage = lecteur.read()["internalMesh"]

        colonnes = {}
        for ax, (nom, x) in zip(axes, LIGNES.items()):
            y_H, u_U0 = profil_vertical(maillage, x)
            ax.plot(u_U0, y_H, color=COULEURS[i], linestyle=STYLES[i],
                    linewidth=2, label=f"itération {t:g}")
            colonnes["y/H"] = y_H
            colonnes[f"u/U0 ({nom})"] = u_U0

        # Sauvegarde des profils bruts de cet instant en CSV
        entete = ",".join(colonnes.keys())
        np.savetxt(f"profils_{t:g}.csv", np.column_stack(list(colonnes.values())),
                   delimiter=",", header=entete, comments="")

    # Mise en forme : axes discrets, grille légère, titres explicites
    for ax, nom in zip(axes, LIGNES):
        ax.set_title(f"Profil vertical à {nom}", fontsize=12)
        ax.set_xlabel("u / U0  (vitesse horizontale adimensionnée)")
        ax.axvline(0, color="#999999", linewidth=0.8)   # u = 0 : séparation aller / retour
        ax.grid(True, color="#e5e5e5", linewidth=0.6)
        ax.set_ylim(0, 1)
        for cote in ("top", "right"):
            ax.spines[cote].set_visible(False)
    axes[0].set_ylabel("y / H  (hauteur adimensionnée)")
    axes[0].legend(frameon=False, loc="best")

    fig.suptitle("Pièce de Nielsen, cas A (k-ε) : convergence des profils de vitesse",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig("profils_convergence.png", dpi=150)
    print("Graphique écrit : profils_convergence.png")

    # Écart maximal entre les deux derniers instants : le chiffre à mettre dans un rapport
    if len(instants) >= 2:
        ecarts = []
        for t_a, t_b in zip(instants[:-1], instants[1:]):
            lecteur.set_active_time_value(t_a)
            m_a = lecteur.read()["internalMesh"]
            lecteur.set_active_time_value(t_b)
            m_b = lecteur.read()["internalMesh"]
            for nom, x in LIGNES.items():
                _, u_a = profil_vertical(m_a, x)
                _, u_b = profil_vertical(m_b, x)
                ecarts.append((t_a, t_b, nom, np.max(np.abs(u_b - u_a))))
        print("\nÉcart maximal de u/U0 entre deux instants successifs :")
        for t_a, t_b, nom, e in ecarts:
            print(f"  {t_a:g} -> {t_b:g}  ({nom}) : {e:.2e}")


if __name__ == "__main__":
    main()
