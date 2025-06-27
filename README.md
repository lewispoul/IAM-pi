# IAM Molecule Viewer

**Une interface moderne pour visualiser et analyser des molécules, faire des calculs de chimie computationnelle, accessible sur Linux et Windows (via WSL).**

---

## 👀 Aperçu rapide

* Visualisation 3D de molécules (ex : methane\_test.xyz fourni)
* Calculs XTB (énergie, HOMO/LUMO, etc.)
* Génération de fichiers .json/.cube (XTB)
* Interface web simple et claire (3Dmol.js)
* Extraction automatique, historique, batch
* Modules avancés : prédiction VoD, DeltaH, extraction PDF, etc.
* 100% open source

---

## 🚦 Prérequis système

### 1. **Linux** (Ubuntu recommandé) OU **Windows 10/11 avec WSL2**

* [Guide d’installation officiel WSL2 (Windows Subsystem for Linux)](https://learn.microsoft.com/fr-fr/windows/wsl/install)
  Ou en un seul copier-coller dans PowerShell (admin) :

  ```powershell
  wsl --install
  ```

  (redémarrer ensuite, choisir Ubuntu)

### 2. **Python 3.8 ou plus**

* Ubuntu : déjà présent ou installer avec :

  ```bash
  sudo apt update && sudo apt install python3 python3-venv python3-pip git
  ```

### 3. **Git**

* Installer avec :

  ```bash
  sudo apt install git
  ```

### 4. **Conda (recommandé)**

* [Miniconda Download – Linux](https://docs.conda.io/en/latest/miniconda.html)
* Installer avec :

  ```bash
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  # Suivre les instructions à l’écran puis redémarrer le terminal
  ```

---

### 📦 Installation à partir d’une archive .zip (alternative)

Si vous avez téléchargé IAM sous forme d’archive `.zip` depuis GitHub, **décompressez-la dans un dossier sans espaces ni accents** (ex : `C:/Users/votre_nom/IAM` ou `~/IAM`), puis ouvrez ce dossier dans votre terminal avant de suivre les instructions ci-dessous.

---

## 🚀 Installation en 5 étapes (la méthode facile)

### **1. Ouvrir un terminal dans Ubuntu (ou WSL/Ubuntu sous Windows)**

### **2. Copier-coller ces commandes** :

```bash
# Télécharger le projet IAM
git clone https://github.com/lewispoul/IAM.git
cd IAM

# Installer Miniconda (si pas déjà fait)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# (Redémarrer le terminal si c'est la première installation conda)
```

### **3. Créer et activer l’environnement de travail**

```bash
conda create -n chem-env python=3.10
conda activate chem-env
```

### **4. Installer toutes les dépendances**

```bash
conda install -c conda-forge rdkit flask flask-cors xtb
```

### **5. Lancer IAM**

```bash
python backend.py
```

Puis ouvrir dans votre navigateur (Windows ou Linux) :
[http://localhost:5000](http://localhost:5000)

---

### 🚀 Option alternative – Installation avec pip/venv

Pour les utilisateurs qui préfèrent ne pas utiliser Conda, vous pouvez utiliser un environnement virtuel Python classique :

```bash
python -m venv iam-env
source iam-env/bin/activate   # Linux/Mac
# ou
iam-env\Scripts\activate      # Windows

pip install -r requirements.txt
```

---

## 🧪 Tester IAM avec un exemple

* Dans l’interface web, cliquez sur “Parcourir”, choisissez `methane_test.xyz` (déjà fourni dans le dossier IAM), puis cliquez sur “Lancer les calculs IAM”.
* Vous verrez un résumé de la molécule et les résultats XTB affichés.

---

## 🗂️ Modules et dossiers principaux

| Dossier                | Rôle principal                                     |
| ---------------------- | -------------------------------------------------- |
| `IAM_Molecule_Engine/` | Calculs moléculaires (XTB, Psi4, etc.)             |
| `IAM_GUI/`             | Interface utilisateur Flask + Web visualisations   |
| `IAM_Knowledge/`       | Datasets, fiches molécules, extraction automatique |
| `IAM_VoD_Predictor/`   | Prédiction performance énergétique (ML/KJ/DFT)     |
| `IAM_Utils/`           | Scripts batch, parsing, fusion, logs, conversions  |
| `IAM_Results/`         | Résultats sauvegardés automatiquement              |

*(Adaptez cette liste selon les dossiers présents dans votre dépôt)*

---

## ❓ Problèmes fréquents et solutions

| Problème                              | Solution                                                                                                                                         |
| ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| **XTB not found**                     | Vérifiez que vous avez bien fait `conda install xtb`                                                                                             |
| **Page web inaccessible**             | Vérifiez que `python backend.py` est bien lancé, ouvrez [http://localhost:5000](http://localhost:5000) dans le navigateur (sur Windows ou Linux) |
| **Permission denied**                 | Ajoutez `sudo` devant la commande si besoin ou vérifiez les droits du dossier                                                                    |
| **Command not found: conda**          | Redémarrez le terminal après installation Miniconda                                                                                              |
| **Cannot open display/3D viewer**     | Essayez un autre navigateur (Chrome ou Edge)                                                                                                     |
| **Erreur sur le port (déjà utilisé)** | Changez le port : `python backend.py --port 5001`                                                                                                |

---

## 📁 Structure du projet (pour comprendre)

```
IAM/
│
├── backend.py               # Démarre l’interface web
├── iam_molecule_engine.py   # Génère, vérifie et convertit les fichiers moléculaires
├── xtb_wrapper.py           # Lancement XTB direct via Python
├── index.html / script.js / style.css    # Interface web moderne
├── methane_test.xyz         # Exemple de molécule pour test
├── reset_iam_env.sh         # Script de nettoyage (optionnel)
└── ...
```

---

## 📝 Astuces et conseils

* **Gardez toujours “chem-env” activé** avant d’utiliser IAM (`conda activate chem-env`)
* **Pour toute question,** contactez [Louis-Philippe Poulin](mailto:votreadresse@exemple.com) ou ouvrez une issue sur GitHub.
* **Vous êtes perdu ?** Refaites chaque étape en copiant-collant les commandes.

---

## 📝 Besoin d’aide ?

En cas de souci, **donnez si possible** :

* 📍 Le nom du script ou module concerné
* 💣 La molécule testée
* 🧪 Le type d’analyse souhaité (XTB, VoD, etc.)
* ⚠️ Le message d’erreur affiché

---

## 📝 À propos du .gitignore et GitHub

Les fichiers temporaires, résultats, environnements virtuels sont automatiquement exclus du dépôt grâce au `.gitignore`. Les actions GitHub automatisées vérifient la bonne configuration avant chaque push.

---

## 🏁 Pour les utilisateurs avancés

* Pour des workflows personnalisés : explorez les fichiers Python (ajout de modules, calculs avancés, intégration ML…).
* IAM peut tourner sur Raspberry Pi, serveurs Linux ou via SSH distant (avec port forwarding).

---

**IAM Molecule Viewer — Québec, Juin 2025**
*Votre assistant pour la chimie computationnelle, prêt à l’emploi !*

---

