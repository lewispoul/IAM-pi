import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from flask import Flask, render_template, request
from IAM_Molecule_Engine.iam_molecule_engine import full_molecule_workflow

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    results = None
    if request.method == 'POST':
        smiles = request.form.get('smiles', '')
        job_name = request.form.get('job_name', 'job')
        if smiles:
            try:
                results = full_molecule_workflow(smiles, job_name)
            except Exception as e:
                results = {'error': str(e)}
    return render_template('index.html', results=results)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)
