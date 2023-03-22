# Get R, Rstudio, python, in ubuntu image based on rocker-version2 
FROM bioconductor/bioconductor_docker:RELEASE_3_16

ENV PATH $PATH:/opt/conda/bin
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc

COPY scraps_conda.yml .
RUN /opt/conda/bin/conda config --add channels defaults && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    conda install -c conda-forge mamba && \
    mamba env update -n base -f scraps_conda.yml && \
    conda clean -afy

RUN R -e 'BiocManager::install("rnabioco/scrapR")'
WORKDIR /home/rstudio
