FROM --platform=amd64 mambaorg/micromamba:1.5.8-noble

COPY environment.yml .
COPY maxquant.zip .

# install and setup of micromamba
USER root
ENV HOME=/root
ENV ENV_NAME=maxquant

RUN echo 'show_banner: false' > ~/.mambarc

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes

RUN  apt-get update -y \
    && apt-get install -y unzip \
    && apt-get clean \
    && unzip maxquant.zip -d /opt/ \
    && ln -s /opt/MaxQuant* /opt/MaxQuant \
    && rm -rf maxquant.zip

# First is necessary for base_image to activate the conda environment second is entrypoint
# which adds the python file to PATH
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
