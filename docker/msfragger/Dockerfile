FROM --platform=amd64 mambaorg/micromamba:2.1.1-ubuntu22.04

WORKDIR /home/mambauser
# Copy environment.yml and requirements.txt
COPY --chown=mambauser:mambauser environment.yml .
COPY --chown=mambauser:mambauser requirements.txt .

# copy over MSFragger bin
COPY --chown=mambauser:mambauser ./MSFragger-4.2.zip /tmp/

# setup of micromamba
USER mambauser
ENV HOME=/home/mambauser
ENV ENV_NAME=medbioinforub-msfragger

RUN echo 'show_banner: false' > ~/.mambarc

RUN micromamba env create -y -f environment.yml \
    && micromamba clean --all --yes

# install / unzip MSFragger
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN unzip /tmp/MSFragger-4.2.zip -d /home/mambauser \
    && rm /tmp/MSFragger-4.2.zip

# First is necessary for base_image to actvate the conda environment second is entrypoint
# which adds the python file to PATH
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
