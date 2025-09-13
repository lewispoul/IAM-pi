"""Flask API étendue pour l'agent ChatGPT IAM avec fonctionnalités IDE."""
import os
import yaml
from flask import Flask, request, jsonify, render_template, Response
from flask_cors import CORS
from flask_socketio import SocketIO, emit
from typing import Optional, Dict, Any
import sys
import subprocess
import threading
import queue
import pty
import select
import termios
import struct
import fcntl
import time

# Add IAM directory to Python path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from IAM_ChatGPT_Integration import IAMChatGPTAgent

app = Flask(__name__, template_folder='templates')
CORS(app)
socketio = SocketIO(app, cors_allowed_origins="*")

# Load config
with open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "agent_config.yaml"), "r", encoding="utf-8") as f:
    config = yaml.safe_load(f)

# Initialize agent
agent = IAMChatGPTAgent(
    api_key=config["gpt4_api_key"],
    model=config.get("gpt4_model", "gpt-4o")
)

# Enable GOD mode
agent.enable_god_mode()

# Store conversations and terminals
conversations: Dict[str, Any] = {}
terminals: Dict[str, Any] = {}

def read_terminal_output(fd, terminal_queue):
    """Read output from terminal and put it in queue."""
    while True:
        try:
            data = os.read(fd, 1024).decode()
            if not data:
                break
            terminal_queue.put(data)
        except:
            break

@app.route('/')
def index():
    """Page principale de l'interface ChatGPT."""
    return render_template('chatgpt.html')

@app.route('/api/chat', methods=['POST'])
def chat():
    """Endpoint API pour le chat."""
    data = request.json
    if not data or 'message' not in data:
        return jsonify({'error': 'Message requis'}), 400
    
    conversation_id = data.get('conversation_id')
    user_message = data['message']
    
    # Get existing conversation or None for new one
    conversation = conversations.get(conversation_id) if conversation_id else None
    
    try:
        # Get response from agent
        response, updated_conversation = agent.chat(user_message, conversation)
        
        # Generate new conversation ID if needed
        if not conversation_id:
            conversation_id = str(len(conversations))
        
        # Store updated conversation
        conversations[conversation_id] = updated_conversation
        
        return jsonify({
            'response': response,
            'conversation_id': conversation_id
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/reset', methods=['POST'])
def reset_conversation():
    """Reset une conversation."""
    data = request.json
    if not data or 'conversation_id' not in data:
        return jsonify({'error': 'Conversation ID requis'}), 400
    
    conversation_id = data['conversation_id']
    if conversation_id in conversations:
        del conversations[conversation_id]
    
    return jsonify({'status': 'success'})

def get_available_port(start_port=5000):
    """Trouve un port disponible."""
    import socket
    port = start_port
    while port < start_port + 100:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(('', port))
                return port
        except OSError:
            port += 1
    raise RuntimeError("Aucun port disponible trouvé")

if __name__ == '__main__':
    # Set environment variables for GOD mode
    os.environ['IAM_GOD_MODE'] = 'TRUE'
    os.environ['IAM_UNLIMITED_ACCESS'] = 'TRUE'
    os.environ['IAM_SYSTEM_COMMANDS'] = 'TRUE'
    
    port = get_available_port()
    print(f"🌐 Interface ChatGPT disponible sur: http://localhost:{port}")
    print("⚡ GOD MODE activé")
    app.run(host='0.0.0.0', port=port, debug=True)