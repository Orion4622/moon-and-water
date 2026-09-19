# La Lune et l'évacuation du moment cinétique

Michel Debailleul, 2026. Licence CC BY-NC 4.0.

Article publié sur Zenodo : https://doi.org/10.5281/zenodo.22844012

## Contenu

| Fichier | Rôle |
|---|---|
| `Lune_article_EN.pdf`, `Lune_article_EN.tex` | Article, version anglaise |
| `Lune_article_FR.pdf`, `Lune_article_FR.tex` | Article, version française |
| `fig_*.png` | Figures (fr et en) |
| `stripping.py` | Moteur d'intégration |
| `analyse.py` | Classement des éléments et statistiques |
| `reproduire.py` | Reproduit tous les nombres et toutes les figures |
| `resultats.json`, `histogramme_jz.npz` | Résultats du calcul |

## Reproduire les calculs

Python 3 avec numpy, scipy et matplotlib :

    python3 reproduire.py              # calcul complet puis figures
    python3 reproduire.py --figures    # figures seules
