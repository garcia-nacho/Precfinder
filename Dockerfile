FROM garcianacho/fhibaseillumina:22022022
LABEL maintainer="Nacho Garcia <iggl@fhi.no>"

COPY Scripts/ /home/docker/Scripts/
COPY Binaries/ /home/docker/Binaries/
RUN chmod +x /home/docker/Binaries/* \
    && chmod +x /home/docker/Scripts/* \
    && mv /home/docker/Binaries/* /usr/bin/
USER docker
WORKDIR /home/docker/Scripts