FROM debian:stretch

RUN apt-get update && \
    apt-get install -y build-essential

RUN apt-get install -y cmake \
                       libgsl0-dev

COPY . /laplus
WORKDIR /laplus
RUN make