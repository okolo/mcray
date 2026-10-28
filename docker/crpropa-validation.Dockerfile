# The local base is built by bin/test-docker.sh from the current mcray sources.
FROM mcray-test-linux-amd64:local
COPY reference/ /reference/
COPY adapters/ /reference-adapters/
RUN cmake -S /reference-adapters -B /reference-build \
      -DCRPROPA_SOURCE_DIR=/reference -DCMAKE_BUILD_TYPE=Release \
    && cmake --build /reference-build --target crpropa_decay_probe mcray_decay_probe -j 4
WORKDIR /src/bin
ENTRYPOINT []
