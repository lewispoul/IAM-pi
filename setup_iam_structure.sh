#!/bin/bash

echo "📁 Création des dossiers IAM..."
mkdir -p ~/IAM/{notebooks,results,data,logs,install,Molecule_Engine/IAM_config,Molecule_Engine/results}

echo "📂 Organisation des fichiers existants..."

# Déplacer les notebooks mal placés
if [ -f ~/IAM/Molecule_Engine/demo.ipynb ]; then
  mv ~/IAM/Molecule_Engine/demo.ipynb ~/IAM/notebooks/
fi

if [ -f ~/IAM/notebooks/iam_molecule_engine.ipynb ]; then
  mv ~/IAM/notebooks/iam_molecule_engine.ipynb ~/IAM/notebooks/_OLD_iam_molecule_engine_backup.ipynb
fi

# Télécharger le .py si absent
if [ ! -f ~/IAM/Molecule_Engine/iam_molecule_engine.py ]; then
  echo "⬇️ Téléchargement du module iam_molecule_engine.py..."
  wget -O ~/IAM/Molecule_Engine/iam_molecule_engine.py "https://raw.githubusercontent.com/openai/sample-projects/main/iam_molecule_engine.py"
fi

# Créer README si absent
if [ ! -f ~/IAM/README.md ]; then
  echo "# IAM – Intelligence Assistée Moléculaire" > ~/IAM/README.md
  echo "Ce projet contient des modules de chimie computationnelle, des notebooks et des outils d'assistance IA pour la modélisation moléculaire." >> ~/IAM/README.md
fi

echo "✅ Structure IAM prête."
