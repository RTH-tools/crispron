FROM tensorflow/tensorflow:latest-gpu

RUN mkdir -p ~/miniconda3 \
    && apt install wget \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh \
    && bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 \
    && rm -rf ~/miniconda3/miniconda.sh \
    && ~/miniconda3/bin/conda init bash \
    && ~/miniconda3/bin/conda init zsh


RUN ~/miniconda3/bin/conda create -y -c bioconda -c conda-forge -n crispron python=3.8
RUN ~/miniconda3/bin/conda install -y -c bioconda -c conda-forge -n crispron biopython viennarna=2.2.5

# activate the env
RUN echo "conda activate crispron" >> ~/.bashrc

# note we might want to mount local data drive
RUN mkdir crispron
WORKDIR /crispron

USER root

LABEL name={NAME}
LABEL version={VERSION}
