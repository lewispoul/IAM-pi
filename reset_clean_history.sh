#!/usr/bin/env bash
set -euo pipefail

# 0) Safety check
git rev-parse --is-inside-work-tree >/dev/null 2>&1 || { echo "Run this inside your repo"; exit 1; }

# 1) Install git-filter-repo if missing
if ! command -v git-filter-repo >/dev/null 2>&1; then
  python3 -m pip install --user git-filter-repo
  export PATH="$HOME/.local/bin:$PATH"
fi

# 2) Make sure we’re on main and up to date (optional)
git fetch origin
git checkout main || git checkout -b main

# 3) Stop tracking secret files in the working tree
echo -e ".env\n*.env\nagent_config.yaml\nstart_iam_with_api.sh\nconfigure_openai.sh" >> .gitignore
git rm -r --cached --ignore-unmatch .env *.env agent_config.yaml start_iam_with_api.sh configure_openai.sh
git add .gitignore
git commit -m "Stop tracking secrets; add .gitignore" || true

# 4) Rewrite history to remove those files everywhere
git filter-repo --force \
  --path .env --path-glob '*.env' \
  --path agent_config.yaml \
  --path start_iam_with_api.sh \
  --path configure_openai.sh \
  --invert-paths

# 5) Push the rewritten history
git remote -v
echo "Force-pushing cleaned history to origin/main…"
git push origin main --force --prune

# 6) Quick verification helper
echo "If you want to verify nothing remains, run:"
echo "  git rev-list --objects --all | grep -E '\\.env|agent_config\\.yaml|start_iam_with_api\\.sh|configure_openai\\.sh' || echo 'No matches ✅'"
