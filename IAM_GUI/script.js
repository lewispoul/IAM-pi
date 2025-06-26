document.addEventListener("DOMContentLoaded", function () {
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
        const viewer = $3Dmol.createViewer("viewer", { backgroundColor: "white" });
        viewer.addModel(contents, "xyz");
        viewer.setStyle({}, { stick: {} });
        viewer.zoomTo();
        viewer.render();
    }

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

    window.launchIAM = async function () {
        const fileInput = document.getElementById('xyzFile');
        if (!fileInput.files.length) {
            alert('Veuillez sélectionner un fichier .xyz');
            return;
        }

        const formData = new FormData();
        formData.append('file', fileInput.files[0]);

        try {
            const response = await fetch('/run_xtb', {
                method: 'POST',
                body: formData
            });

            const result = await response.json();

            document.getElementById('summaryContent').textContent = JSON.stringify({
                "Total Energy (Eh)": result["total energy"],
                "HOMO-LUMO gap (eV)": result["HOMO-LUMO gap/eV"],
                "Dipole moment (D)": result["dipole"],
                "XTB version": result["xtb version"]
            }, null, 2);

            document.getElementById('inputFileContent').textContent = result["input_file"] || "Aucun input fourni.";
            document.getElementById('outputFileContent').textContent = result["output_file"] || JSON.stringify(result, null, 2);

            showTab('summary');
        } catch (err) {
            alert('Erreur réseau ou de calcul : ' + err);
        }
    };
});
