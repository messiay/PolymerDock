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
        "Run this notebook step-by-step to execute the *real* physics simulation (GNINA docking + OpenMM 10ns MD) on a free Colab GPU.\n",
        "\n",
        "**Important:** Run each cell in order. Cell 1 will restart the kernel automatically — that's normal!"
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
        "print('✅ Conda environment is active.')\n",
        "\n",
        "# Determine the CONDA Python path (sys.executable on Colab is the system Python, NOT conda)\n",
        "import os\n",
        "CONDA_PYTHON = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'python')\n",
        "print(f'Conda Python: {CONDA_PYTHON}')\n",
        "print(f'sys.executable: {__import__(\"sys\").executable}  ← this is the SYSTEM Python, we do NOT use this')\n",
        "print()\n",
        "\n",
        "!conda install -y -c conda-forge openmm mdtraj rdkit biopython openmmforcefields openff-toolkit ambertools\n",
        "print('\\n✅ Conda packages installed.')"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 3. Install pip packages INTO the CONDA Python (not the system one!)\n",
        "import os\n",
        "CONDA_PYTHON = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'python')\n",
        "print(f'Installing pip packages via: {CONDA_PYTHON}')\n",
        "\n",
        "!{CONDA_PYTHON} -m pip install --quiet meeko gemmi pdbfixer fastapi uvicorn pydantic matplotlib pyyaml pandas scipy numpy\n",
        "\n",
        "# Download docking engines & Cloudflare tunnel binary\n",
        "!wget -q https://github.com/gnina/gnina/releases/download/v1.0.3/gnina -O /usr/bin/gnina\n",
        "!chmod +x /usr/bin/gnina\n",
        "!wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/bin/vina\n",
        "!chmod +x /usr/bin/vina\n",
        "!wget -q https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-linux-amd64 -O /usr/bin/cloudflared\n",
        "!chmod +x /usr/bin/cloudflared\n",
        "print('\\n✅ Pip packages and system binaries installed.')"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 4. Verify OpenMM + GPU environment BEFORE launching the server\n",
        "# We test using the CONDA Python, which is what the server will use.\n",
        "import os, subprocess\n",
        "CONDA_PYTHON = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'python')\n",
        "\n",
        "verify_script = '''\n",
        "import sys, importlib\n",
        "print(f\"Python executable: {sys.executable}\")\n",
        "print(f\"Python version:    {sys.version}\")\n",
        "print()\n",
        "\n",
        "packages = [\"openmm\", \"mdtraj\", \"rdkit\", \"openmmforcefields\", \"openff.toolkit\"]\n",
        "all_ok = True\n",
        "for pkg in packages:\n",
        "    try:\n",
        "        mod = importlib.import_module(pkg)\n",
        "        ver = getattr(mod, \"__version__\", getattr(mod, \"version\", \"?\"))\n",
        "        print(f\"  ✅ {pkg:25s} → {ver}\")\n",
        "    except ImportError as e:\n",
        "        print(f\"  ❌ {pkg:25s} → MISSING: {e}\")\n",
        "        all_ok = False\n",
        "\n",
        "print()\n",
        "try:\n",
        "    import openmm\n",
        "    for i in range(openmm.Platform.getNumPlatforms()):\n",
        "        p = openmm.Platform.getPlatform(i)\n",
        "        print(f\"  Platform {i}: {p.getName()}\")\n",
        "    try:\n",
        "        cuda = openmm.Platform.getPlatformByName(\"CUDA\")\n",
        "        print(f\"\\\\n  🚀 CUDA platform available! MD will run on GPU.\")\n",
        "    except Exception:\n",
        "        print(f\"\\\\n  ⚠️  CUDA not available, will fall back to CPU (slower but still real MD).\")\n",
        "except Exception as e:\n",
        "    print(f\"  ❌ OpenMM platform check failed: {e}\")\n",
        "\n",
        "if all_ok:\n",
        "    print(\"\\\\n🎉 All packages verified — REAL MD simulations will run!\")\n",
        "else:\n",
        "    print(\"\\\\n⚠️  Some packages are missing. Re-run Cell 2.\")\n",
        "'''\n",
        "\n",
        "# Run the verification inside the CONDA Python\n",
        "result = subprocess.run([CONDA_PYTHON, '-c', verify_script], capture_output=True, text=True)\n",
        "print(result.stdout)\n",
        "if result.stderr:\n",
        "    print('STDERR:', result.stderr)"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 5. Retrieve project source files\n",
        "# Clones from your GitHub repository directly. If cloning fails (e.g. private repo restrictions),\n",
        "# it falls back to unzipping 'SimDock_Project.zip' if you uploaded it manually.\n",
        "import os\n",
        "if not os.path.exists('polymerDock'):\n",
        "    !git clone https://github.com/messiay/PolymerDock.git polymerDock || (unzip -q -o SimDock_Project.zip -d polymerDock && echo 'Fallback to local SimDock_Project.zip successful')\n",
        "else:\n",
        "    print('polymerDock directory already exists, skipping clone.')\n",
        "%cd polymerDock"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 6. Launch the SimDock Web Application with Cloudflare Tunnel!\n",
        "import subprocess, time, re, os\n",
        "\n",
        "# Use the CONDA Python — NOT sys.executable!\n",
        "CONDA_PYTHON = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'python')\n",
        "\n",
        "# Build subprocess environment with conda paths\n",
        "env = os.environ.copy()\n",
        "conda_prefix = os.environ['CONDA_PREFIX']\n",
        "conda_bin = os.path.join(conda_prefix, 'bin')\n",
        "conda_lib = os.path.join(conda_prefix, 'lib')\n",
        "env['PATH'] = conda_bin + ':' + env.get('PATH', '')\n",
        "env['LD_LIBRARY_PATH'] = conda_lib + ':' + env.get('LD_LIBRARY_PATH', '')\n",
        "\n",
        "print(f'Conda Python: {CONDA_PYTHON}')\n",
        "print(f'CONDA_PREFIX: {conda_prefix}')\n",
        "print()\n",
        "\n",
        "print('Starting FastAPI Uvicorn Server...')\n",
        "server_log = open('server.log', 'w')\n",
        "proc = subprocess.Popen(\n",
        "    [CONDA_PYTHON, '-m', 'uvicorn', 'api:app', '--host', '0.0.0.0', '--port', '8501'],\n",
        "    env=env,\n",
        "    stdout=server_log,\n",
        "    stderr=subprocess.STDOUT\n",
        ")\n",
        "time.sleep(5)\n",
        "\n",
        "# Check if server started successfully\n",
        "if proc.poll() is not None:\n",
        "    server_log.close()\n",
        "    print('❌ Server crashed on startup! Logs:')\n",
        "    with open('server.log') as f:\n",
        "        print(f.read())\n",
        "else:\n",
        "    print('✅ Server is running on port 8501.')\n",
        "\n",
        "print('\\nStarting Cloudflare Tunnel...')\n",
        "cf_log_path = 'cloudflare.log'\n",
        "with open(cf_log_path, 'w') as f:\n",
        "    subprocess.Popen(['cloudflared', 'tunnel', '--url', 'http://localhost:8501'], stdout=f, stderr=f)\n",
        "\n",
        "# Poll the log file for the public tunnel URL\n",
        "tunnel_url = None\n",
        "print('Retrieving public HTTPS URL from Cloudflare...')\n",
        "for _ in range(15):\n",
        "    time.sleep(1)\n",
        "    if os.path.exists(cf_log_path):\n",
        "        with open(cf_log_path, 'r') as f:\n",
        "            content = f.read()\n",
        "            match = re.search(r'https://[a-zA-Z0-9-]+\\.trycloudflare\\.com', content)\n",
        "            if match:\n",
        "                tunnel_url = match.group(0)\n",
        "                break\n",
        "\n",
        "if tunnel_url:\n",
        "    print('\\n✅ SUCCESS! Click the link below to open your app in a new tab:')\n",
        "    print(f'\\n👉 {tunnel_url} 👈\\n')\n",
        "else:\n",
        "    print('\\n⚠️ Could not automatically detect Cloudflare URL.')\n",
        "    print('Check cloudflare.log for details, or try the port window fallback:')\n",
        "    try:\n",
        "        from google.colab import output\n",
        "        output.serve_kernel_port_as_window(8501)\n",
        "    except Exception:\n",
        "        pass"
      ],
      "outputs": [],
      "execution_count": None
    },
    {
      "cell_type": "code",
      "metadata": {},
      "source": [
        "# 7. (Debug) View server logs if something goes wrong\n",
        "# Run this cell if the dashboard says 'OpenMM unavailable' or the server crashed.\n",
        "with open('server.log') as f:\n",
        "    print(f.read())"
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
