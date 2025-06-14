#!/bin/bash

echo "🧼 Nettoyage de l'environnement IAM..."

# Supprimer les anciens résultats
rm -rf ~/IAM/results/*
rm -f ~/IAM/*.log ~/IAM/*.out

# Supprimer les fichiers cube si existants
find ~/IAM -name "*.cube" -delete

# Réinitialiser Jupyter Kernel
jupyter kernelspec remove python3 -f
source ~/IAM/venv/bin/activate
python3 -m ipykernel install --user --name=python3

echo "✅ Nettoyage terminé. Lancement du serveur Jupyter..."
jupyter lab
