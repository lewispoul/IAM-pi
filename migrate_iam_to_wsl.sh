#!/bin/bash

# 📁 Dossier source (actuel)
SRC_DIR="$(pwd)"

# 📁 Dossier destination sous WSL
DEST_DIR="$HOME/IAM"

echo "🚀 Migration de IAM vers $DEST_DIR"

# Créer le dossier destination s’il n’existe pas
mkdir -p "$DEST_DIR"

# Copier tous les fichiers et sous-dossiers sauf .git, __pycache__, .vscode, .DS_Store
rsync -av --progress "$SRC_DIR/" "$DEST_DIR/" \
    --exclude='.git' \
    --exclude='__pycache__' \
    --exclude='.vscode' \
    --exclude='.DS_Store'

echo "✅ Migration terminée."
