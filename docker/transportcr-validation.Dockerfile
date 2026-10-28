FROM mcray-test-linux-amd64:local
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    libxerces-c-dev python3 git && rm -rf /var/lib/apt/lists/*
COPY transportcr/ /reference/transportcr/
ENV CMAKE_BUILD_PARALLEL_LEVEL=2
# The stock script ends its build invocation with a bare --parallel, which
# overrides CMAKE_BUILD_PARALLEL_LEVEL. Keep the script and cap only this call.
RUN bash -c 'cmake() { if [[ "$1" == "--build" && "${@: -1}" == "--parallel" ]]; then command cmake "$@" 2; else command cmake "$@"; fi; }; export -f cmake; bash /reference/transportcr/bin/build.sh Release' 
WORKDIR /src/bin
ENTRYPOINT ["bash", "/src/tests/transportcr/run_psb_comparison.sh", "/reference/transportcr", "/results"]
