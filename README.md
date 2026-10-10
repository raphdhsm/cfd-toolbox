# cfd-toolbox

Cas de validation et outils CFD (OpenFOAM) de **Raphaël Duhesme**, ingénieur en mécanique des fluides
(INSA Hauts-de-France, spécialité mécanique et énergétique).

Chaque cas suit la même démarche : un benchmark avec données de référence, un maillage paramétré en
Python, une validation chiffrée contre les mesures, une étude de convergence en maillage.

## Outils

- [OpenFOAM 13](https://openfoam.org) : solveur
- [classy_blocks](https://github.com/damogranlabs/classy_blocks) : maillages structurés paramétrés en Python
- [foamlib](https://github.com/gerlero/foamlib) : pilotage des cas en Python (clonage, paramètres, lancement)
- [pyvista](https://github.com/pyvista/pyvista) : post-traitement automatisé
- ParaView : visualisation

## Cas

| Cas | Domaine | Résultat clé |
| --- | --- | --- |
| [Pièce ventilée de Nielsen (IEA Annex 20)](cas/piece_nielsen_IEA_annex20) | Ventilation / HVAC | k-ε : écart moyen 3 à 7 % de U₀ aux mesures ; k-ω SST 1,7 à 2,9 fois moins précis ; GCI 2 à 3 % |
