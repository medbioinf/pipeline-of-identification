# AMD64 needed explicitly on ARM as some software is only available in AMD64
FROM mambaorg/micromamba:2.3.1-ubuntu24.04

# Copy backend and environment.yml
COPY --chown=mambauser:mambauser environment.yml /home/mambauser/
COPY --chown=mambauser:mambauser requirements.txt /home/mambauser/

# setup build environment (currently needed for psm-utils forks)
USER root
RUN apt update \
    && apt install -y libglib2.0-0 build-essential cmake \
    && apt clean

# setup of micromamba
USER mambauser
ENV HOME=/home/mambauser
ENV ENV_NAME=medbioinforub-ident-comparison
WORKDIR /home/mambauser

RUN echo 'show_banner: false' > ~/.mambarc \
    && micromamba env create -y -f /home/mambauser/environment.yml \
    && micromamba clean --all --yes

# entrypoint for micromamba
ENTRYPOINT [ "/usr/local/bin/_entrypoint.sh"]
