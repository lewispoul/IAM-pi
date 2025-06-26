#!/bin/bash

echo "🔧 Installation des dépendances système (OpenBabel, swig)..."
sudo apt update
sudo apt install -y swig openbabel libopenbabel-dev

echo "📦 Activation de l'environnement virtuel..."
source ~/IAM/venv/bin/activate

echo "📦 Installation de openbabel-wheel dans le venv..."
pip install openbabel-wheel

echo "📦 Installation de pybel (wrapper Python pour Open Babel)..."
pip install openbabel

echo "✅ Vérification de l'import Python..."
python -c "import pybel; print('✔️ pybel importé depuis :', pybel.__file__)"
