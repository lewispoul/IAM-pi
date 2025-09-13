# Endpoints Flask pour calculator
# Import du module
from GeneratedScripts.calculator import *

from flask import Flask, request, jsonify

app = Flask(__name__)

@app.route('/calculate', methods=['POST'])
def calculate():
    data = request.json
    expression = data.get('expression', '')
    
    try:
        # Évaluer l'expression mathématique
        result = str(eval(expression))
        return jsonify({'result': result}), 200
    except Exception as e:
        return jsonify({'error': 'Invalid Input'}), 400

@app.route('/clear', methods=['POST'])
def clear():
    # Réinitialiser l'expression
    return jsonify({'result': ''}), 200

if __name__ == '__main__':
    app.run(debug=True)