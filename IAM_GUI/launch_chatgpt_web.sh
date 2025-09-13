#!/bin/bash

# Couleurs pour le terminal
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${YELLOW}🚀 Lancement de l'interface web ChatGPT IAM${NC}"
echo "=================================================="

# Vérification de l'environnement conda
if ! command -v conda &> /dev/null; then
    echo -e "${RED}❌ Conda n'est pas installé${NC}"
    exit 1
fi

# Activation de l'environnement conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate chem-env || conda activate nox-env

# Vérification des dépendances
echo -e "\n${YELLOW}🔍 Vérification des dépendances...${NC}"
python3 -c "import flask" 2>/dev/null || pip install flask
python3 -c "import flask_cors" 2>/dev/null || pip install flask-cors
python3 -c "import flask_socketio" 2>/dev/null || pip install flask-socketio
python3 -c "import eventlet" 2>/dev/null || pip install eventlet

# Variables d'environnement pour le GOD MODE
export IAM_GOD_MODE=TRUE
export IAM_UNLIMITED_ACCESS=TRUE
export IAM_SYSTEM_COMMANDS=TRUE

# Lancement du serveur Flask
cd "$(dirname "$0")"
echo -e "\n${GREEN}🌐 Démarrage du serveur web...${NC}"
python3 chatgpt_api.py