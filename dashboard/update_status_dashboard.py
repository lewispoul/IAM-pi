from datetime import datetime

status_data = {
    "XTB": "✅ Terminé",
    "Psi4": "🔧 En cours",
    "VoD (Kamlet–Jacobs)": "✅ Terminé",
    "VoD (ML)": "🔧 En cours",
    "Keshavarz Benchmark": "🕒 À faire",
    "IAM Interface (UI)": "✅ Terminé",
    "Cube Viewer": "🔧 En cours",
    "IAM Agent": "✅ Terminé",
    "Copilot Memory": "🕒 À faire",
    "Dashboard": "✅ Terminé"
}

output_path = "IAM_StatusDashboard.html"

html = f"""<!DOCTYPE html>
<html lang="fr">
<head>
    <meta charset="UTF-8">
    <title>IAM Dashboard</title>
    <style>
        body {{
            font-family: sans-serif;
            background-color: #f2f2f2;
            padding: 20px;
        }}
        h1 {{
            color: #333;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            background-color: white;
        }}
        th, td {{
            padding: 10px;
            border: 1px solid #ccc;
            text-align: left;
        }}
        .status-ok {{ color: green; font-weight: bold; }}
        .status-pending {{ color: orange; font-weight: bold; }}
        .status-todo {{ color: gray; font-weight: bold; }}
    </style>
</head>
<body>
    <h1>📊 Tableau de bord IAM</h1>
    <p>🕓 Dernière mise à jour : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    <table>
        <tr><th>Module</th><th>Statut</th></tr>
"""

for module, status in status_data.items():
    cls = "status-ok" if "✅" in status else "status-pending" if "🔧" in status else "status-todo"
    html += f"<tr><td>{module}</td><td class='{cls}'>{status}</td></tr>\n"

html += """
    </table>
</body>
</html>
"""

with open(output_path, "w", encoding="utf-8") as f:
    f.write(html)
print(f"✅ Dashboard mis à jour : {output_path}")
