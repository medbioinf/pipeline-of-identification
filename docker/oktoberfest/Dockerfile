# AMD64 needed explicitly on ARM as some sofwtare is only available in AMD64
FROM --platform=amd64 mambaorg/micromamba:2.1.0-ubuntu22.04

WORKDIR /home/mambauser
# Copy backend and environment.yml
COPY --chown=mambauser:mambauser environment.yml .
COPY --chown=mambauser:mambauser requirements.txt .

USER root

RUN apt update \
    && apt install -y libglib2.0-0 git \
    && apt clean

USER mambauser
ENV HOME=/home/mambauser
ENV ENV_NAME=oktoberfest

RUN echo 'show_banner: false' > ~/.mambarc

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes

# TODO: remove build-essential

USER root
# First is necessary for base_image to actvate the conda environment second is entrypoint
# which adds the python file to PATH
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
