FROM --platform=amd64 mambaorg/micromamba:1.5.8-noble

WORKDIR /home/mambauser

COPY --chown=mambauser:mambauser environment.yml .
COPY ./maxquant.zip ./maxquant.zip

# install and setup of micromamba and openms
USER mambauser
ENV HOME=/home/mambauser
ENV ENV_NAME=maxquant

RUN echo 'show_banner: false' > ~/.mambarc

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes

USER root

RUN  apt-get update -y \
    && apt-get install -y unzip \
    && apt-get clean \
    && unzip maxquant.zip -d /opt/ \
    && ln -s /opt/MaxQuant* /opt/MaxQuant \
    && rm -rf maxquant.zip

# First is necessary for base_image to activate the conda environment second is entrypoint
# which adds the python file to PATH
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
