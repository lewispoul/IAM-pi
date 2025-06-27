document.addEventListener("DOMContentLoaded", function () {
    console.log('DOMContentLoaded fired, DOM is ready.');
    function showTab(tabName) {
        document.querySelectorAll('.tab-content').forEach(tab => {
            tab.classList.remove('active');
        });
        document.querySelectorAll('.tab-button').forEach(btn => {
            btn.classList.remove('active');
        });

        document.getElementById(tabName).classList.add('active');
        document.querySelector(`[onclick="showTab('${tabName}')"]`).classList.add('active');
    }

    window.showTab = showTab;

    function renderMolecule(contents) {
        const viewerDiv = document.getElementById("viewer");
        if (!viewerDiv) {
            alert('3D viewer element (id="viewer") not found in DOM.');
            return;
        }
        let viewer;
        try {
            viewer = $3Dmol.createViewer(viewerDiv, { backgroundColor: "white" });
        } catch (e) {
            alert('3Dmol.js failed to initialize: ' + e);
            return;
        }
        viewer.addModel(contents, "xyz");
        viewer.setStyle({}, { stick: {} });
        viewer.zoomTo();
        viewer.render();
    }

    window.renderMolecule = renderMolecule;

    document.getElementById('xyzFile').addEventListener('change', function (event) {
        const file = event.target.files[0];
        if (file && file.name.endsWith('.xyz')) {
            const reader = new FileReader();
            reader.onload = function (e) {
                renderMolecule(e.target.result);
            };
            reader.readAsText(file);
        }
    });

    // Show/hide Psi4 options based on method
    window.togglePsi4Options = function () {
        const method = document.getElementById('method').value;
        document.getElementById('psi4Options').style.display = (method === 'psi4') ? '' : 'none';
    };

    // Load molecule from pasted .xyz/.mol content
    window.loadFromPaste = function () {
        const contents = document.getElementById('xyzPaste').value.trim();
        if (!contents) {
            alert('Please paste .xyz or .mol content.');
            return;
        }
        renderMolecule(contents);
    };

    window.launchIAM = async function () {
        alert('launchIAM called');
        const fileInput = document.getElementById('xyzFile');
        const pasteInput = document.getElementById('xyzPaste').value.trim();
        console.log('fileInput.files:', fileInput.files);
        if (!fileInput.files.length && !pasteInput) {
            alert('No file selected and no pasted content.');
        }
        if (fileInput.files.length) {
            alert('File selected: ' + fileInput.files[0].name);
        }
        if (pasteInput) {
            alert('Using pasted content.');
        }
        let formData = new FormData();
        if (fileInput.files.length) {
            formData.append('file', fileInput.files[0]);
        } else if (pasteInput) {
            // Convert pasted text to a Blob and append as file
            const blob = new Blob([pasteInput], { type: 'text/plain' });
            formData.append('file', blob, 'pasted.xyz');
        } else {
            alert('Veuillez sélectionner un fichier .xyz ou coller le contenu dans la zone de texte.');
            return;
        }
        formData.append('method', document.getElementById('method').value);
        formData.append('basis', document.getElementById('basis').value);
        formData.append('charge', document.getElementById('charge').value);
        formData.append('multiplicity', document.getElementById('multiplicity').value);
        formData.append('calcType', document.getElementById('calcType').value);
        formData.append('solvent', document.getElementById('solvent').value);
        if (document.getElementById('method').value === 'psi4') {
            formData.append('functional', document.getElementById('functional').value);
        }
        try {
            const response = await fetch('/run_xtb', {
                method: 'POST',
                body: formData
            });
            const result = await response.json();

            // Format summary output for XTB and future Psi4
            if (result.success && result.xtb_json) {
                const xtb = result.xtb_json;
                let html = '<table class="summary-table">';
                if (xtb["total energy"]) html += `<tr><td>Total Energy (Eh)</td><td>${xtb["total energy"]}</td></tr>`;
                if (xtb["HOMO-LUMO gap/eV"]) html += `<tr><td>HOMO-LUMO gap (eV)</td><td>${xtb["HOMO-LUMO gap/eV"]}</td></tr>`;
                if (xtb["dipole"]) html += `<tr><td>Dipole moment (D)</td><td>${xtb["dipole"]}</td></tr>`;
                if (xtb["xtb version"]) html += `<tr><td>XTB version</td><td>${xtb["xtb version"]}</td></tr>`;
                html += '</table>';
                // Collapsible section for arrays
                if (xtb["charges"]) {
                    html += `<details><summary>Atomic Charges</summary><pre>${JSON.stringify(xtb["charges"], null, 2)}</pre></details>`;
                }
                if (xtb["wbo"]) {
                    html += `<details><summary>Wiberg Bond Orders</summary><pre>${JSON.stringify(xtb["wbo"], null, 2)}</pre></details>`;
                }
                document.getElementById('summaryContent').innerHTML = html;
            } else if (result.success && result.psi4_json) {
                const psi4 = result.psi4_json;
                let html = '<table class="summary-table">';
                if (psi4.final_energy !== undefined) html += `<tr><td>Final Energy (Eh)</td><td>${psi4.final_energy}</td></tr>`;
                if (psi4.psi4_version) html += `<tr><td>Psi4 Version</td><td>${psi4.psi4_version}</td></tr>`;
                html += '</table>';
                html += `<details><summary>Stdout</summary><pre>${psi4.stdout || ''}</pre></details>`;
                html += `<details><summary>Stderr</summary><pre>${psi4.stderr || ''}</pre></details>`;
                document.getElementById('summaryContent').innerHTML = html;
            } else if (result.details) {
                document.getElementById('summaryContent').textContent = result.details;
            } else {
                document.getElementById('summaryContent').textContent = JSON.stringify(result, null, 2);
            }

            // Always clear input file content before updating
            document.getElementById('inputFileContent').textContent = '';
            if (result.file_preview) {
                document.getElementById('inputFileContent').textContent = result.file_preview;
            } else if (result.details && result.details.includes('Aucun fichier reçu')) {
                document.getElementById('inputFileContent').textContent = 'No file received by backend.';
            } else {
                document.getElementById('inputFileContent').textContent = 'No file';
            }
            document.getElementById('outputFileContent').textContent = result["output_file"] || JSON.stringify(result, null, 2);

            showTab('summary');
        } catch (err) {
            alert('Erreur réseau ou de calcul : ' + err);
        }
    };

    // Load molecule from SMILES input and render in 3D viewer
    window.loadFromSMILES = async function () {
        const smiles = document.getElementById('smilesInput').value.trim();
        if (!smiles) {
            alert('Please enter a SMILES string.');
            return;
        }
        // Call backend to convert SMILES to XYZ
        const response = await fetch('/smiles_to_xyz', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles })
        });
        const data = await response.json();
        if (data.success && data.xyz) {
            renderMolecule(data.xyz);
        } else {
            alert('Failed to convert SMILES to 3D structure.');
        }
    };

    // Load molecule from Ketcher (2D editor)
    window.loadFromKetcher = async function () {
        const ketcherFrame = document.getElementById('ketcherFrame');
        const ketcher = ketcherFrame.contentWindow;
        if (!ketcher) {
            alert('Ketcher not loaded.');
            return;
        }
        // Get molfile from Ketcher
        ketcher.postMessage({ type: 'get-molfile' }, '*');
        window.addEventListener('message', async function handler(event) {
            if (event.data && event.data.type === 'molfile') {
                window.removeEventListener('message', handler);
                const molfile = event.data.molfile;
                // Call backend to convert molfile to XYZ
                const response = await fetch('/molfile_to_xyz', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ molfile })
                });
                const data = await response.json();
                if (data.success && data.xyz) {
                    renderMolecule(data.xyz);
                } else {
                    alert('Failed to convert 2D structure to 3D.');
                }
            }
        });
    };

    // Placeholder for searchMolecule
    window.searchMolecule = function () {
        alert('Search functionality not implemented yet.');
    };
});
