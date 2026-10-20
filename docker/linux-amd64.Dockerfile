ARG UBUNTU_VERSION=24.04
FROM ubuntu:${UBUNTU_VERSION}
RUN test "$(dpkg --print-architecture)" = amd64 \
    && apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
       ca-certificates cmake make gcc g++ gfortran libgsl-dev \
       libboost-filesystem-dev libboost-system-dev libyaml-cpp-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY src/ src/
COPY tests/ tests/
COPY examples/ examples/
COPY bin/build.sh bin/build.sh
COPY bin/tables/ bin/tables/
RUN bash bin/build.sh Release
WORKDIR /src/bin
ENTRYPOINT ["ctest", "--test-dir", "/src/build/Release", "--output-on-failure"]
