FROM mambaorg/micromamba:2.3.1-ubuntu24.04

# Copy environment.yml and requirements.txt
COPY --chown=mambauser:mambauser environment.yml /home/mambauser/
COPY --chown=mambauser:mambauser requirements.txt /home/mambauser/

# copy over MSFragger bin
COPY --chown=mambauser:mambauser ./MSFragger-4.2.zip /tmp/

# setup of micromamba
USER mambauser
ENV HOME=/home/mambauser
ENV ENV_NAME=medbioinforub-msfragger
WORKDIR /home/mambauser

RUN echo 'show_banner: false' > ~/.mambarc \
    && micromamba env create -y -f /home/mambauser/environment.yml \
    && micromamba clean --all --yes

# install / unzip MSFragger
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN unzip /tmp/MSFragger-4.2.zip -d /home/mambauser \
    && rm /tmp/MSFragger-4.2.zip

# First is necessary for base_image to actvate the conda environment second is entrypoint
# which adds the python file to PATH
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
