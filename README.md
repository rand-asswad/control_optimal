# Contrôle Optimal

Ce répertoire est un mini-projet de contrôle optimal
dans le cadre de ma formation de mathématiques appliquées.

Dans ce mini-projet, un problème de contrôle est traité
mathématiquement dans le rapport (HTML et PDF),
et une implémentation `python` est fournit à l'aide
des libraires `sympy`, `numpy` et `matplotlib`.

# Installation

Télécharger/cloner le répertoire et naviguer au répertoire.

```sh
git clone https://github.com/rand-asswad/controle_optimal.git
cd controle_optimal
```

Afin de lancer le code localement, il est recommandé
d'installer les paquets python nécessaire dans
un environnement virtuel.

```sh
# Créer un environnement virtuel
python3 -m venv venv

# Activer l'environnement virtuel
source venv/bin/activate

# Installer les paquets requis
pip install -r requirements.txt
```

Lancer le programme
```sh
python3 main.py
```

# Rapport

## R Markdown

Le rapport est généré à l'aide de
[R Markdown](https://rmarkdown.rstudio.com/)
et hébergé sur [Github Pages](https://pages.github.com/).

R Markdown est un outil génial de 
- [knitr](https://yihui.org/knitr/): une librairie **R**
  pour exécuter des parties du code dans les fichiers Markdown.
- [pandoc](https://pandoc.org/): un paquet qui permet
  de convertir des fichiers d'un format à un autre.

Pour plus de détails, consulter la documentation:
- [R Markdown: The Definitive Guide](https://bookdown.org/yihui/rmarkdown/)
- [bookdown](https://bookdown.org/yihui/bookdown/)

## Installation

Installer les paquets **R** localement sur la console `r`

```r
install.packages('bookdown')
install.packages('tinytex')
tinytex::install_tinytex()
```

## Génération du rapport

```sh
# Générer book.pdf (via LaTeX)
make pdf

# Générer index.html
make html
```