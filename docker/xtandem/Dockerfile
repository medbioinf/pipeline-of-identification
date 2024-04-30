FROM --platform=amd64 ubuntu:24.04

ENV DEBIAN_FRONTEND noninteractive
ENV COMET_VERSION=${COMET_VERSION}

RUN apt-get update -y \
    && apt-get install -y curl \
    && apt-get clean

COPY docker/entrypoint.sh /usr/local/bin

# install X!Tandem
RUN apt-get update -y \
    && apt-get install -y unzip \
    && apt-get clean

RUN curl -L -o /tmp/tandem.zip ftp://ftp.thegpm.org/projects/tandem/source/tandem-linux-17-02-01-4.zip \
    && unzip /tmp/tandem.zip tandem-linux-17-02-01-4/bin/static_link_ubuntu/tandem.exe -d /tmp/ \
    && mv /tmp/tandem-linux-17-02-01-4/bin/static_link_ubuntu/tandem.exe /usr/local/bin/tandem.exe  \
    && chmod a+x /usr/local/bin/tandem.exe \
    && ln -sf /usr/local/bin/tandem.exe /usr/local/bin/tandem \
    && rm -rf /tmp/tandem.zip /tmp/tandem-linux-17-02-01-4

ENTRYPOINT [ "/usr/local/bin/entrypoint.sh" ]