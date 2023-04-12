FROM eddelbuettel/r2u:focal

# install R packages
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        time \
        libxml2-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libgsl-dev \
        libhdf5-dev && \
    install.r Seurat tidyverse patchwork ggrastr rjson gridExtra && \
    install.r BiocManager remotes && \
    installGithub.r mojaveazure/seurat-disk

#RUN ln -s /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/local/bin/ && \
#    installBioc.r BSgenome.Hsapiens.UCSC.hg38
RUN ln -s /usr/local/lib/R/site-library/littler/examples/installBioc.r /usr/bin/ && \
    installBioc.r BSgenome.Hsapiens.UCSC.hg38


RUN install.r devtools hdf5r RSpectra && \
    MAKEFLAGS="-j 4" r -e 'devtools::install_github(c("bnprks/BPCells")); quit(status=!requireNamespace("BPCells"))'


# install miniconda 
RUN apt-get update -q && \
    apt-get install -q -y --no-install-recommends \
        bzip2 \
        ca-certificates \
        git \
        libglib2.0-0 \
        libsm6 \
        libxext6 \
        libxrender1 \
        mercurial \
        openssh-client \
        procps \
        subversion \
        wget \
        zstd \
        gawk \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENV PATH /opt/conda/bin:$PATH

CMD [ "/bin/bash" ]

# Leave these args here to better use the Docker build cache
ARG CONDA_VERSION=py310_22.11.1-1

RUN set -x && \
    UNAME_M="$(uname -m)" && \
    if [ "${UNAME_M}" = "x86_64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-x86_64.sh"; \
        SHA256SUM="00938c3534750a0e4069499baf8f4e6dc1c2e471c86a59caa0dd03f4a9269db6"; \
    elif [ "${UNAME_M}" = "s390x" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-s390x.sh"; \
        SHA256SUM="a150511e7fd19d07b770f278fb5dd2df4bc24a8f55f06d6274774f209a36c766"; \
    elif [ "${UNAME_M}" = "aarch64" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-aarch64.sh"; \
        SHA256SUM="48a96df9ff56f7421b6dd7f9f71d548023847ba918c3826059918c08326c2017"; \
    elif [ "${UNAME_M}" = "ppc64le" ]; then \
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-${CONDA_VERSION}-Linux-ppc64le.sh"; \
        SHA256SUM="4c86c3383bb27b44f7059336c3a46c34922df42824577b93eadecefbf7423836"; \
    fi && \
    wget "${MINICONDA_URL}" -O miniconda.sh -q && \
    echo "${SHA256SUM} miniconda.sh" > shasum && \
    if [ "${CONDA_VERSION}" != "latest" ]; then sha256sum --check --status shasum; fi && \
    mkdir -p /opt && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh shasum && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    find /opt/conda/ -follow -type f -name '*.a' -delete && \
    find /opt/conda/ -follow -type f -name '*.js.map' -delete && \
    /opt/conda/bin/conda clean -afy

# ensure at least 4GB RAM when building the container
RUN conda install -c dranew bcl2fastq
RUN conda install -c conda-forge poppler 
RUN conda install -c bioconda -c conda-forge snakemake
RUN conda install -c bioconda -c conda-forge bowtie2 fastp subread star umi_tools pysam samtools htslib 

# fix a samtools dependency issue by reinstalling ncurses
RUN conda install -c conda-forge ncurses

RUN install.r Cairo && \
    r -e 'BiocManager::install("TFBSTools", type="source")' && \
    installGithub.r GreenleafLab/chromVAR GreenleafLab/chromVARmotifs && \
    r -e 'devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories(), \
          dependency=TRUE, type="source")' && \
    r -e 'library(ArchR); ArchR::installExtraPackages()'
