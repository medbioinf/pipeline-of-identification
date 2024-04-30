FROM --platform=amd64 ubuntu:24.04

ENV DEBIAN_FRONTEND noninteractive
ENV COMET_VERSION=${COMET_VERSION}

RUN apt-get update -y \
    && apt-get install -y curl \
    && apt-get clean

COPY docker/entrypoint.sh /usr/local/bin

# install Percolator
ARG PERCOLATOR_REL_VERSION=rel-3-06-05
ARG PERCOLATOR_SHORT_VERSION=3-06

RUN apt-get update -y \
    && apt-get install -y libgomp1 \
    && apt-get clean

RUN curl -L -o /tmp/percolator.deb https://github.com/percolator/percolator/releases/download/${PERCOLATOR_REL_VERSION}/percolator-v${PERCOLATOR_SHORT_VERSION}-linux-amd64.deb \
    && dpkg -i /tmp/percolator.deb \
    && apt-get install -f -y \
    && rm -rf /tmp/percolator.deb \
    && apt-get clean

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]