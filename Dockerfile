# Get R, Rstudio, python, in ubuntu image based on rocker-version2 
FROM bioconductor/bioconductor_docker:RELEASE_3_16

ENV PATH $PATH:/opt/conda/bin
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

COPY scraps_conda.yml .
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    conda install -c conda-forge mamba && \
    mamba env update -n base -f scraps_conda.yml && \
    conda clean -afy

# fail the build if the solve did not produce Snakemake 9 (scraps needs >= 9, Python >= 3.11)
RUN snakemake --version && python -c "import sys; assert sys.version_info >= (3, 11), sys.version"

RUN R -e 'BiocManager::install("rnabioco/scrapR")'
WORKDIR /home/rstudio
