from flask import Flask, request, jsonify, render_template
from flask_cors import CORS
import os
import subprocess
import tempfile
import json

app = Flask(__name__)
CORS(app)

@app.route('/')
def index():
    return render_template('iam_viewer_connected.html')

@app.route('/run_xtb', methods=['POST'])
def run_xtb():
    if "file" not in request.files:
        return jsonify({"success": False, "details": "Aucun fichier reçu"}), 400

    xyz_file = request.files["file"]
    if xyz_file.filename == "":
        return jsonify({"success": False, "details": "Nom de fichier vide"}), 400

    try:
        with tempfile.TemporaryDirectory() as tempdir:
            xyz_path = os.path.join(tempdir, "molecule.xyz")
            xyz_file.save(xyz_path)

            xtb_command = ["xtb", xyz_path, "--opt", "--json", "--gfn", "2"]
            result = subprocess.run(xtb_command, cwd=tempdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            json_path = os.path.join(tempdir, "xtbout.json")
            if not os.path.exists(json_path):
                return jsonify({"success": False, "details": "Fichier xtbout.json non trouvé"}), 500

            with open(json_path, "r") as f:
                xtb_data = json.load(f)

            return jsonify({"success": True, "xtb_json": xtb_data})

    except Exception as e:
        return jsonify({"success": False, "details": str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000)