FROM eclipse-temurin:23-jre-noble

# copy over MSFragger bin
COPY ./MSFragger-4.1.zip /tmp/

# install dependencies
RUN apt update \
    && apt install -y unzip \
    && rm -rf /var/lib/apt/lists/*

# install MS-GF+ and bash
RUN unzip /tmp/MSFragger-4.1.zip -d /opt/ \
    && rm /tmp/MSFragger-4.1.zip

CMD ["java", "-jar", "/opt/MSFragger-4.1/MSFragger-4.1.jar"]
