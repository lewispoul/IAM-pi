# ⚗️ README IAM – Installation et Premier Lancement

## 📦 Décompression

1. Télécharge et **extrais l’archive .zip** dans le dossier : `C:\Users\pouli\OneDrive\Bureau\IAM` *(ou tout autre dossier sans espace ni accent)*

---

## 🐍 Création de l’environnement Python

### Option 1 – Conda (recommandé)

```bash
conda env create -f chem-env.yaml
conda activate chem-env
```

*Si tu n’as pas conda, télécharge Miniconda ou Anaconda.*

### Option 2 – Pip (si tu préfères)

```bash
python -m venv iam-env
iam-env\Scripts\activate       # Windows
# ou
source iam-env/bin/activate   # Linux/Mac

pip install -r requirements.txt
```

---

## ⚙️ Installation IAM

1. (Optionnel) Pour faciliter la vie en terminal :

```bash
chmod +x iam.sh
# (optionnel) ajouter à ton PATH
sudo ln -s $(pwd)/iam.sh /usr/local/bin/iam
```

1. **Initialise la base** :

```bash
python iam_update_db.py
# OU via le menu interactif :
./iam.sh menu
```

---

## 🌐 Lancer l’interface web IAM

```bash
cd IAM_GUI
python app.py
# Puis ouvre http://127.0.0.1:5000 dans ton navigateur (Chrome/Firefox)
```

Le mode sombre est automatique, ou activable dans l’UI.

---

## 🖥️ Menu CLI

```bash
./iam.sh menu
```

Ou :

```bash
iam menu
```

---

## 📂 Dossiers principaux

* `IAM_Molecule_Engine/` — core calcul XTB/Psi4
* `IAM_GUI/` — interface web moderne (Ketcher, historique, presets…)
* `IAM_VoD_Predictor/` — prédiction performance EM
* `IAM_Knowledge/` — datasets, fiches molécules, extraction auto, rapports
* `IAM_Utils/` — scripts batch, extraction, fusion, logs
* `IAM_Results/` — résultats calculs, historiques

---

## 💡 Premiers tests

* Lance un calcul sur le nitrométhane (preset DFT ou XTB)
* Essaie l’historique interactif, l’ajout batch, l’extraction de datasets
* Teste le bouton de mise à jour auto de la base
* Joue avec IAM-Copilot dès qu’il sera activé dans l’UI

---

## ⚠️ Si tu rencontres le moindre souci

* Envoie-moi ici : le message d’erreur, l’étape concernée, le script/dossier impliqué
* Je t’envoie un patch ou le correctif sur-le-champ (aucun risque de rester bloqué !)

---

## 🚀 Prêt pour la suite

* Dès le téléchargement, tu peux commencer à jouer/tester tout ce que tu veux !
* Tu peux enrichir la base (ajout de PDF, CSV…), IAM fait tout automatiquement

---

**Définition – **\`\`** :**

Un fichier `.gitignore` sert à **exclure certains fichiers ou dossiers** du suivi par Git. Typiquement, tu y mets :

* les fichiers temporaires (`*.pyc`, `__pycache__/`)
* les environnements (`env/`, `.venv/`)
* les résultats de calcul ou fichiers volumineux (`*.log`, `IAM_Results/`)
* les fichiers de configuration personnels (`.vscode/`, `config.yaml`)

Cela évite d’encombrer ton dépôt avec des fichiers inutiles ou sensibles.

---

**Prêt à être copié directement dans ton GitHub.**
