#!/bin/bash

echo "🔧 Setting up IAM project environment..."

# 1. Create IAM project folder
mkdir -p ~/IAM/{core/modules,canvas,chat_history,data,logs,tools}
cd ~/IAM

# 2. Initialize Git
git init
touch README.md
echo "IAM – AI Chemistry Assistant (initialized on $(date))" > README.md

# 3. Create virtual environment
sudo apt update
sudo apt install -y python3-venv python3-pip
python3 -m venv venv
source venv/bin/activate

# 4. Install base Python packages
pip install --upgrade pip
pip install numpy pandas openbabel rdkit

# 5. Install Open Babel CLI
sudo apt install -y openbabel

# 6. Download xtb binary (optional manual step if unavailable in repo)
# echo "📦 xtb must be manually downloaded and built if not available on Pi"

# 7. Create starter runner.py
cat <<EOL > core/runner.py
# IAM Runner – Core execution script (placeholder)
print("IAM Runner launched")
EOL

# 8. Create .gitignore
cat <<EOF > .gitignore
venv/
__pycache__/
*.log
*.tmp
*.pyc
EOF

echo "✅ IAM setup complete. Run 'source ~/IAM/venv/bin/activate' to begin."
