# La Lune et l'évacuation du moment cinétique

Michel Debailleul, 2026. Licence CC BY-NC 4.0.

## Contenu

| Fichier | Rôle |
|---|---|
| `Lune_article_FR.tex`, `Lune_article_FR.pdf` | Article, version française |
| `Lune_article_EN.tex`, `Lune_article_EN.pdf` | Article, version anglaise |
| `fig_capture_*.png`, `fig_controle_*.png`, `fig_bilan_*.png` | Figures (`fr` et `en`) |
| `stripping.py` | Moteur d'intégration : trajectoire de B, critère de soulèvement, particules-tests |
| `analyse.py` | Classement des éléments et statistiques |
| `reproduire.py` | Reproduit tous les nombres et toutes les figures de l'article |
| `resultats.json`, `histogramme_jz.npz` | Résultats produits par `reproduire.py` |

## Compiler l'article

Les figures doivent se trouver dans le même dossier que le `.tex`.

    pdflatex Lune_article_FR.tex
    pdflatex Lune_article_FR.tex

## Reproduire les calculs

Python 3 avec numpy, scipy et matplotlib.

    python3 reproduire.py              # calcul complet puis figures (quelques minutes)
    python3 reproduire.py --figures    # figures seules, à partir de resultats.json
