"""
Maillage de la pièce ventilée de Nielsen (IEA Annex 20, cas 2D) avec classy_blocks.

Ce script écrit le fichier system/blockMeshDict d'un cas OpenFOAM.
Il sait générer deux maillages, à partir de la même géométrie :

    A : pour le k-epsilon avec lois de paroi (première maille ~2,5 cm aux parois)
    B : pour le k-omega SST sans loi de paroi (première maille ~1,5 mm, y+ ~ 1)

Utilisation (dans WSL) :
    python3 maillage_nielsen.py A cas_kepsilon
    python3 maillage_nielsen.py B cas_komegaSST
    python3 maillage_nielsen.py A cas_fin 1.5      # 3e argument optionnel : facteur de raffinement

Puis, dans le dossier du cas :
    blockMesh
    checkMesh
"""

import sys
from pathlib import Path

import classy_blocks as cb

# =============================================================================
# 1. GÉOMÉTRIE (en mètres) — valeurs de la spécification de Nielsen (1990)
# =============================================================================
H = 3.0       # hauteur de la pièce
L = 3.0 * H   # longueur de la pièce (L/H = 3)
h = 0.056 * H # hauteur de la fente d'entrée (0,168 m)
t = 0.16 * H  # hauteur de la sortie (0,48 m)
EP = 0.1      # épaisseur dans la direction z (cas 2D : une seule maille, valeur sans importance)

# Hauteurs des coupures entre les blocs (voir le schéma des 3 bandes)
Y_COUPURE_BAS = t        # changement sortie -> mur sur le mur droit
Y_COUPURE_HAUT = H - h   # changement mur -> entrée sur le mur gauche (= couche de cisaillement du jet)

# =============================================================================
# 2. PARAMÈTRES DE MAILLAGE — c'est ici que les maillages A et B diffèrent
# =============================================================================
# Toutes les tailles sont en mètres. Une "taille" = hauteur (ou largeur) d'une maille.
#
#   paroi       : taille de la première maille contre les murs (pilote le y+)
#   cisaillement: taille des mailles de part et d'autre de la coupure du haut (y = H - h),
#                 là où le jet rencontre l'air immobile
#   coeur       : taille maximale des mailles au centre de la pièce
#   expansion   : rapport de taille maximal entre deux mailles voisines (règle : <= 1,2)
PARAMETRES = {
    "A": {"paroi": 0.025, "cisaillement": 0.008, "coeur": 0.08, "expansion": 1.15},
    "B": {"paroi": 0.0015, "cisaillement": 0.004, "coeur": 0.08, "expansion": 1.15},
}


def decouper_avec_progression(bloc, axe, longueur, taille_debut, taille_fin, taille_max, expansion):
    """Découpe un côté de bloc en mailles qui grandissent depuis chaque extrémité.

    On part de 'taille_debut' à une extrémité et de 'taille_fin' à l'autre.
    Les mailles grandissent avec le taux 'expansion' (progression géométrique)
    jusqu'à 'taille_max', puis restent constantes au milieu si la place le permet.

    classy_blocks accepte plusieurs "chops" sur le même côté : chacun couvre une
    fraction de la longueur (length_ratio), ce qui permet une progression en plusieurs zones.
    """
    r = expansion

    def longueur_croissance(a, s):
        # Longueur occupée par des mailles qui passent de la taille a à la taille s
        # avec un rapport r entre voisines (somme d'une suite géométrique)
        return max(0.0, (s * r - a) / (r - 1)) if s > a else 0.0

    l_debut = longueur_croissance(taille_debut, taille_max)
    l_fin = longueur_croissance(taille_fin, taille_max)

    if l_debut + l_fin < longueur:
        # Cas 1 : assez de place pour une zone uniforme au milieu
        zones = [
            (l_debut, taille_debut, taille_max),
            (longueur - l_debut - l_fin, taille_max, taille_max),
            (l_fin, taille_max, taille_fin),
        ]
    else:
        # Cas 2 : côté trop court pour atteindre taille_max (ex. la fente).
        # On cherche la taille s où les deux progressions se rejoignent.
        s = (longueur * (r - 1) + taille_debut + taille_fin) / (2 * r)
        if s <= max(taille_debut, taille_fin):
            # Pas de "bosse" possible : simple progression d'un bout à l'autre
            zones = [(longueur, taille_debut, taille_fin)]
        else:
            l_d = longueur_croissance(taille_debut, s)
            zones = [(l_d, taille_debut, s), (longueur - l_d, s, taille_fin)]

    for longueur_zone, t_debut, t_fin in zones:
        if longueur_zone > 1e-9:
            bloc.chop(axe, length_ratio=longueur_zone / longueur, start_size=t_debut, end_size=t_fin)


def construire_maillage(version, facteur=1.0):
    """Construit le maillage A ou B.

    facteur : coefficient de raffinement pour l'étude de convergence en maillage.
      facteur > 1 raffine (1,5 = mailles 1,5 fois plus petites), facteur < 1 déraffine.
      La taille de la PREMIÈRE MAILLE À LA PAROI n'est pas modifiée : avec des lois de paroi,
      la changer ferait varier y+ et mélangerait deux effets (maillage et modèle de paroi).
    """
    p = dict(PARAMETRES[version])
    p["coeur"] /= facteur
    p["cisaillement"] /= facteur
    p["expansion"] = p["expansion"] ** (1.0 / facteur)   # progression plus douce si on raffine
    mesh = cb.Mesh()

    # -------------------------------------------------------------------------
    # 3. LES TROIS BLOCS (bandes horizontales)
    # cb.Box(coin_min, coin_max) crée un bloc rectangulaire.
    # Noms des faces dans classy_blocks :
    #   axe x (0) : left / right    axe y (1) : front / back    axe z (2) : bottom / top
    # -------------------------------------------------------------------------
    bloc_bas = cb.Box([0, 0, 0], [L, Y_COUPURE_BAS, EP])
    bloc_milieu = cb.Box([0, Y_COUPURE_BAS, 0], [L, Y_COUPURE_HAUT, EP])
    bloc_haut = cb.Box([0, Y_COUPURE_HAUT, 0], [L, H, EP])

    # -------------------------------------------------------------------------
    # 4. DÉCOUPAGE EN MAILLES
    # -------------------------------------------------------------------------
    # Direction x (axe 0) : commune aux trois blocs (classy_blocks la propage
    # automatiquement aux blocs voisins). Fin près des murs gauche et droit,
    # plus gros au milieu.
    decouper_avec_progression(bloc_bas, 0, L, p["paroi"], p["paroi"], p["coeur"], p["expansion"])

    # Direction y (axe 1), bloc par bloc, de bas en haut :
    # - bloc du bas : fin au sol, puis on remonte vers la taille "cœur" à y = t
    decouper_avec_progression(
        bloc_bas, 1, Y_COUPURE_BAS, p["paroi"], p["coeur"] / 2, p["coeur"], p["expansion"]
    )
    # - bloc du milieu : raccord avec le bloc du bas en bas, fin en haut (couche de cisaillement)
    decouper_avec_progression(
        bloc_milieu, 1, Y_COUPURE_HAUT - Y_COUPURE_BAS,
        p["coeur"] / 2, p["cisaillement"], p["coeur"], p["expansion"],
    )
    # - bloc du haut (la fente) : fin en bas (cisaillement), "paroi" au plafond
    decouper_avec_progression(
        bloc_haut, 1, H - Y_COUPURE_HAUT, p["cisaillement"], p["paroi"], p["coeur"], p["expansion"]
    )

    # Direction z (axe 2) : une seule maille, c'est un cas 2D
    for bloc in (bloc_bas, bloc_milieu, bloc_haut):
        bloc.chop(2, count=1)

    # -------------------------------------------------------------------------
    # 5. CONDITIONS LIMITES (patches)
    # -------------------------------------------------------------------------
    bloc_haut.set_patch("left", "inlet")     # la fente d'entrée, en haut du mur gauche
    bloc_bas.set_patch("right", "outlet")    # la sortie, en bas du mur droit
    for bloc in (bloc_bas, bloc_milieu, bloc_haut):
        bloc.set_patch(["bottom", "top"], "frontAndBack")  # faces avant/arrière du cas 2D

    for bloc in (bloc_bas, bloc_milieu, bloc_haut):
        mesh.add(bloc)

    # Toutes les faces extérieures non nommées (sol, plafond, murs) -> patch "walls"
    mesh.set_default_patch("walls", "wall")
    mesh.modify_patch("inlet", "patch")
    mesh.modify_patch("outlet", "patch")
    mesh.modify_patch("frontAndBack", "empty")
    return mesh


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4) or sys.argv[1] not in PARAMETRES:
        print(__doc__)
        sys.exit(1)

    version, dossier_cas = sys.argv[1], Path(sys.argv[2])
    facteur = float(sys.argv[3]) if len(sys.argv) == 4 else 1.0
    (dossier_cas / "system").mkdir(parents=True, exist_ok=True)

    mesh = construire_maillage(version, facteur)
    sortie = dossier_cas / "system" / "blockMeshDict"
    mesh.write(str(sortie))
    print(f"Maillage {version} écrit dans {sortie}")

    # blockMesh a besoin d'un controlDict pour se lancer. On en écrit un minimal
    # s'il n'existe pas encore (on le complétera pour le calcul à l'étape suivante).
    control = dossier_cas / "system" / "controlDict"
    if not control.exists():
        control.write_text(
            "FoamFile\n{\n    format      ascii;\n    class       dictionary;\n"
            "    object      controlDict;\n}\n\n"
            "application     foamRun;\nsolver          incompressibleFluid;\n\n"
            "startFrom       startTime;\nstartTime       0;\nstopAt          endTime;\n"
            "endTime         2000;\ndeltaT          1;\n"
            "writeControl    timeStep;\nwriteInterval   200;\n"
        )
        print(f"controlDict minimal écrit dans {control}")
