import json
import zipfile
import os
import subprocess
import sys

# Step 0: Build React App if Node/npm is present and frontend/ exists
frontend_path = os.path.join(os.path.dirname(__file__), 'frontend')
if os.path.exists(frontend_path):
    print("Detected frontend/ directory. Building React production bundle...")
    try:
        # Check npm command exists
        import shutil
        if shutil.which("npm") is not None:
            # Run npm run build
            print("Running 'npm run build' inside frontend/...")
            res = subprocess.run("npm run build", shell=True, cwd=frontend_path, check=True)
            if res.returncode == 0:
                print("React production bundle built successfully in frontend/dist.")
        else:
            print("Warning: npm is not available on PATH. Skipping build. Make sure frontend/dist is built manually.")
    except Exception as e:
        print(f"Warning: Failed to run build command: {e}. Packaging existing files.")

# Create Notebook
notebook = {
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {},
      "source": [
        "# 🧬 SimDock Polymer v2.0 - Cloud GPU Environment\n",
        "Run this notebook step-by-step to execute the *real* physics simulation (GNINA docking + OpenMM 10ns MD) on a free Colab GPU."
      ]
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 1. Setup Conda (The kernel will automatically restart after this finishes. That is normal!)\n",
        "!pip install -q condacolab\n",
        "import condacolab\n",
        "condacolab.install()"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 2. Install heavy chemistry dependencies via Conda\n",
        "import condacolab\n",
        "condacolab.check()\n",
        "!conda install -y -c conda-forge openmm mdtraj rdkit biopython openmmforcefields openff-toolkit"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 3. Install Python pip packages & System Binaries\n",
        "!pip install meeko gemmi pdbfixer fastapi uvicorn pydantic matplotlib pyyaml pandas scipy numpy\n",
        "\n",
        "!wget -q https://github.com/gnina/gnina/releases/download/v1.0.3/gnina -O /usr/bin/gnina\n",
        "!chmod +x /usr/bin/gnina\n",
        "!wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/bin/vina\n",
        "!chmod +x /usr/bin/vina\n",
        "!wget -q https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-linux-amd64 -O /usr/bin/cloudflared\n",
        "!chmod +x /usr/bin/cloudflared"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 4. Retrieve project source files\n",
        "# Clones from your GitHub repository directly. If cloning fails (e.g. private repo restrictions),\n",
        "# it falls back to unzipping 'SimDock_Project.zip' if you uploaded it manually.\n",
        "!git clone https://github.com/messiay/PolymerDock.git polymerDock || (unzip -q -o SimDock_Project.zip -d polymerDock && echo 'Fallback to local SimDock_Project.zip successful')\n",
        "%cd polymerDock"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 5. Launch the SimDock Web Application!\n",
        "import sys\n",
        "import subprocess\n",
        "from google.colab import output\n",
        "import time\n",
        "\n",
        "print(\"Starting FastAPI Uvicorn Server in Conda-Colab environment...\")\n",
        "# Launch backend+frontend on port 8501 using sys.executable (points to active Conda-Colab kernel)\n",
        "subprocess.Popen([sys.executable, \"-m\", \"uvicorn\", \"api:app\", \"--host\", \"0.0.0.0\", \"--port\", \"8501\"])\n",
        "time.sleep(3)\n",
        "\n",
        "print(\"\\n✅ SUCCESS! Click the link below to open your app in a new tab:\")\n",
        "output.serve_kernel_port_as_window(8501)"
      ],
      "outputs": [],
      "execution_count": None
    }
  ],
  "metadata": {
    "accelerator": "GPU",
    "colab": {
      "gpuType": "T4"
    },
    "kernelspec": {
      "display_name": "Python 3",
      "name": "python3"
    }
  },
  "nbformat": 4,
  "nbformat_minor": 0
}

with open("SimDock_Colab.ipynb", "w") as f:
    json.dump(notebook, f, indent=2)

print("Created SimDock_Colab.ipynb")

# Create Zip file
def zipdir(path, ziph):
    for root, dirs, files in os.walk(path):
        # Exclude unnecessary/heavy folders and node_modules
        if any(exclude in root for exclude in ['.git', '__pycache__', 'results', '.gemini', 'node_modules']):
            continue
        for file in files:
            # Exclude the zip itself and the notebook
            if file.endswith('.zip') or file.endswith('.ipynb'):
                continue
            
            file_path = os.path.join(root, file)
            arcname = os.path.relpath(file_path, path)
            ziph.write(file_path, arcname)

print("Packaging project files...")
with zipfile.ZipFile('SimDock_Project.zip', 'w', zipfile.ZIP_DEFLATED) as zipf:
    zipdir('.', zipf)

print("Created SimDock_Project.zip")
