# Runtime image for profiling: the web slicer backend (Ubuntu 24.04, xvfb) plus perf.
FROM mirdain_web_slicer-backend:latest
USER root
RUN apt-get update && apt-get install -y --no-install-recommends linux-tools-generic binutils && rm -rf /var/lib/apt/lists/* \
 && ln -sf "$(ls -d /usr/lib/linux-tools/*/ | head -1)perf" /usr/local/bin/perf-any
