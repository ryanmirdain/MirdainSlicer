ARG BASE_IMAGE=ubuntu:24.04
FROM ${BASE_IMAGE}

ARG CMAKE_VERSION=4.3.0

ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y sudo ca-certificates curl tar

COPY build_linux.sh /tmp/orcaslicer/build_linux.sh
COPY scripts/linux.d /tmp/orcaslicer/scripts/linux.d

WORKDIR /tmp/orcaslicer

RUN chmod +x ./build_linux.sh
RUN ./build_linux.sh -ur
RUN if [[ -n "${CMAKE_VERSION}" ]] ; then \
        case "$(uname -m)" in \
            x86_64|amd64) cmake_arch="x86_64" ;; \
            aarch64|arm64) cmake_arch="aarch64" ;; \
            *) cmake_arch="" ;; \
        esac ; \
        if [[ -n "${cmake_arch}" ]] ; then \
            cmake_root="/opt/cmake-${CMAKE_VERSION}-linux-${cmake_arch}" ; \
            if ! curl -fsSL "https://github.com/Kitware/CMake/releases/download/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}-linux-${cmake_arch}.tar.gz" | tar -xz -C /opt ; then \
                echo "Warning: failed to install CMake ${CMAKE_VERSION}; falling back to the distro cmake package." ; \
            elif [[ -d "${cmake_root}" ]] ; then \
                ln -sf "${cmake_root}/bin/"* /usr/local/bin/ ; \
            fi ; \
        else \
            echo "Skipping GitHub Actions CMake install for unsupported architecture $(uname -m)." ; \
        fi ; \
    fi
RUN rm -rf /var/lib/apt/lists/* /tmp/orcaslicer
