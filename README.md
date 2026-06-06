# 🧬 SimDock Polymer v2.0

Universal Polymer-Enzyme Catalytic Simulation Engine mapping target enzymes against polymer chains using docking, coordinate-driven growth, geometry filtering, and Molecular Dynamics.

---

## 🚀 Execution Guide

### 1. Running in Google Colab (Recommended)
You can run this pipeline inside a free Google Colab GPU instance using the provided notebook:
1. Open a new Google Colab notebook and change your runtime type to **GPU** (**Runtime > Change runtime type > GPU**).
2. Upload the `SimDock_Colab.ipynb` file from this project to your Colab session.
3. Run the cells step-by-step.
4. Copy the public IP address displayed in Step 5's output, click the `loca.lt` link, paste the IP, and submit!

> [!WARNING]  
> **Localtunnel Connection Blips:** When serving Streamlit over `localtunnel` in Google Colab, you might occasionally see red boxes in the UI saying:
> `TypeError: Failed to fetch dynamically imported module: ...`
> This is a known localtunnel forwarding limitation where requests for Streamlit's static widget JavaScript bundles are dropped. 
> **Solution:** Simply perform a **hard reload of your browser tab (Ctrl+F5 or Shift+Reload)** to force localtunnel to fetch the missing asset files.

---

## 🛠️ Docker Setup
To run the simulation engine inside a fully self-contained Docker container:

```bash
# Build the Docker image
docker build -t polymerdock .

# Run the container mapping Streamlit's port
docker run -p 8501:8501 polymerdock
```

Open your browser and navigate to `http://localhost:8501`.
