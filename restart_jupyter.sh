#!/bin/bash

echo "[🧹] Nettoyage de Jupyter..."
pkill -9 -f jupyter
rm -rf ~/.local/share/jupyter/runtime

echo "[🚀] Démarrage de JupyterLab dans venv..."
cd ~/IAM
source venv/bin/activate
jupyter lab
