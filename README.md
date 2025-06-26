# ⚗️ IAM – Personal Chemistry AI Assistant

[![.gitignore check](https://github.com/lewispoul/IAM/actions/workflows/check_gitignore.yml/badge.svg)](https://github.com/lewispoul/IAM/actions/workflows/check_gitignore.yml)

---

## 📦 Installation et Premier Lancement

### 1. Décompression

Télécharge et **extrait l’archive `.zip`** dans un dossier sans espace ni accent (ex. `C:\Users\pouli\OneDrive\Bureau\IAM` ou `~/IAM`).

---

### 2. Création de l’environnement Python

#### ✅ Option 1 – Conda (recommandé)

```bash
conda env create -f chem-env.yaml
conda activate chem-env
```

💡 *Installe [Miniconda](https://docs.conda.io/en/latest/miniconda.html) si tu ne l’as pas.*

#### ✅ Option 2 – Pip

```bash
python -m venv iam-env
iam-env\Scripts\activate       # Windows
# ou
source iam-env/bin/activate   # Linux/Mac

pip install -r requirements.txt
```

---

### 3. Initialisation IAM (Terminal)

```bash
chmod +x iam.sh
# (optionnel) pour un accès global :
sudo ln -s $(pwd)/iam.sh /usr/local/bin/iam

# Démarrer
./iam.sh menu
# ou
iam menu
```

---

## 🌐 Lancer l’interface Web IAM

```bash
cd IAM_GUI
python app.py
```

Ensuite ouvre : [http://127.0.0.1:5000](http://127.0.0.1:5000)

💡 Mode sombre automatique.

---

## 🧪 Premiers Tests

- [x] Calcul XTB sur le **nitrométhane**
- [x] Génération des fichiers `.json`, `.cube`, etc.
- [x] Visualisation 3D avec 3Dmol.js
- [x] Prédiction VoD et deltaH
- [x] Test de l’historique et extraction automatique depuis PDF

---

## 🧠 Modules Principaux

| Dossier                   | Rôle principal                                      |
|--------------------------|-----------------------------------------------------|
| `IAM_Molecule_Engine/`   | Calculs moléculaires (XTB, Psi4, etc.)              |
| `IAM_GUI/`               | Interface utilisateur Flask + Web visualisations    |
| `IAM_Knowledge/`         | Datasets, fiches molécules, extraction automatique  |
| `IAM_VoD_Predictor/`     | Prédiction performance énergétique (ML/KJ/DFT)      |
| `IAM_Utils/`             | Scripts batch, parsing, fusion, logs, conversions   |
| `IAM_Results/`           | Résultats sauvegardés automatiquement               |

---

## 🔁 GitHub Synchronisation

1. Ton dépôt est déjà connecté (origin : `https://github.com/lewispoul/IAM.git`)
2. Le fichier `.gitignore` est vérifié automatiquement :
   [![.gitignore check](https://github.com/lewispoul/IAM/actions/workflows/check_gitignore.yml/badge.svg)](https://github.com/lewispoul/IAM/actions/workflows/check_gitignore.yml)

3. Pour pousser les changements :

```bash
git add .
git commit -m "💬 Mise à jour IAM"
git push origin main
```

---

## 🧾 Définition : `.gitignore`

Le `.gitignore` permet d’exclure :

- Les fichiers temporaires Python (`__pycache__/`, `*.pyc`)
- Les environnements (`env/`, `.venv/`)
- Les résultats de calcul (`IAM_Results/`, `*.log`)
- Les dossiers système (`.ipynb_checkpoints/`, `.vscode/`, `.DS_Store`)

💡 Cela évite de "polluer" le dépôt avec des fichiers locaux non pertinents.

---

## 💡 Conseils

- Garde un œil sur les résultats dans `IAM_Results/`
- Utilise le bouton “Mise à jour” pour enrichir la base automatiquement
- Commence par les molécules simples pour valider le pipeline (ex: CH₃NO₂, NH₃, TATB, etc.)

---

## 🤝 Besoin d'aide ?

Tu peux poser une question en spécifiant :

- 📍 Le nom du script
- 💣 La molécule testée (si applicable)
- 🧪 Le type d’analyse souhaité (XTB, VoD, etc.)
- ⚠️ Le message d’erreur s’il y a lieu

---

IAM est un assistant scientifique extensible — **chaque calcul, chaque fiche, chaque script peut être amélioré, enrichi ou automatisé.** 🧬

---
