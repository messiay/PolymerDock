# Use a Miniconda base image
FROM condaforge/miniforge3:latest

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    libglu1-mesa \
    libxi6 \
    libxrender1 \
    libxrandr2 \
    libxcursor1 \
    libxcomposite1 \
    libasound2t64 \
    libdbus-1-3 \
    libfontconfig1 \
    libxkbcommon0 \
    && rm -rf /var/lib/apt/lists/*

# Install GNINA (statically linked Linux binary)
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.0.3/gnina -O /usr/bin/gnina && \
    chmod +x /usr/bin/gnina

# Install AutoDock Vina
RUN wget -q https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64 -O /usr/bin/vina && \
    chmod +x /usr/bin/vina

# Set up the work directory
WORKDIR /app

# Copy requirements
COPY requirements.txt /app/requirements.txt

# Create conda environment with key packages that are best installed via conda
RUN conda install -y -c conda-forge \
    openmm \
    mdtraj \
    rdkit \
    biopython \
    openmmforcefields \
    openff-toolkit \
    && conda clean -afy

# Install the rest of python dependencies
RUN pip install --no-cache-dir \
    meeko \
    gemmi \
    pdbfixer \
    fastapi \
    uvicorn \
    pydantic \
    matplotlib \
    pyyaml \
    pandas \
    scipy \
    numpy

# Copy the rest of the application
COPY . /app

# Expose web app port
EXPOSE 8501

# Run the FastAPI application using Uvicorn
CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8501"]
