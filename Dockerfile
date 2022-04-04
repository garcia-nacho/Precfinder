FROM garcianacho/fhibaseillumina:04042022
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

RUN Rscript -e "install.packages(c('keras', 'entropy', 'abind'))"
USER docker
RUN Rscript -e "keras::install_keras()"

USER root
COPY Scripts/ /home/docker/Scripts/
COPY Binaries/ /home/docker/Binaries/
COPY Models/ /Models/
RUN chmod +x /home/docker/Binaries/* \
    && chmod +x /home/docker/Scripts/* \
    && chmod -R 777 /Models \
    && mv /home/docker/Binaries/* /usr/bin/
USER docker
WORKDIR /home/docker/Scripts
