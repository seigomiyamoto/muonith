# =============================================================================
# MUONITH — Installation verification Dockerfile
# =============================================================================
#
# Builds MUONITH in an isolated Ubuntu container using g++ and OpenBLAS.
# Dependencies (Eigen3, nlohmann-json, fmt, spdlog) are fetched automatically
# via CMake FetchContent — no manual library installation is required.
#
# --- Selecting the Ubuntu version --------------------------------------------
#
#   # Default (Ubuntu 24.04 LTS — recommended):
#   docker build -f docker/Dockerfile -t muonith-test .
#
#   # Ubuntu 22.04 LTS (GCC 13 PPA and Kitware CMake added automatically):
#   docker build -f docker/Dockerfile --build-arg UBUNTU_VERSION=22.04 -t muonith-test .
#
# --- Run interactively -------------------------------------------------------
#
#   docker run -it muonith-test bash
#
#   # Inside the container, the built binaries are in
#   # /workspace/build-debug/ and /workspace/build-release/
#
# --- Clean up ----------------------------------------------------------------
#
#   docker rmi muonith-test
#
# =============================================================================

# Global ARGs — must appear before the first FROM for use in FROM directives.
ARG UBUNTU_VERSION=24.04
ARG GIFSKI_VERSION=1.34.0

# ---------------------------------------------------------------------------
# Stage 1: Build gifski (high-quality GIF encoder)
# Uses a separate stage so the Rust toolchain does not bloat the final image.
# ---------------------------------------------------------------------------
FROM rust:slim AS gifski-builder
ARG GIFSKI_VERSION
RUN cargo install gifski --version ${GIFSKI_VERSION}

# ---------------------------------------------------------------------------
# Stage 2: Main build
# ---------------------------------------------------------------------------
FROM ubuntu:${UBUNTU_VERSION}

ENV DEBIAN_FRONTEND=noninteractive

# Re-declare ARG after FROM so it is available in subsequent RUN commands.
ARG UBUNTU_VERSION

# ---------------------------------------------------------------------------
# System packages
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      ninja-build \
      ccache \
      pkg-config \
      git \
      time \
      python3 \
      python3-pip \
      software-properties-common \
      wget \
      gpg \
      gpg-agent \
      ca-certificates \
      libopenblas-openmp-dev \
      liblapacke-dev \
      libtbb-dev \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# gifski binary (copied from builder stage)
# ---------------------------------------------------------------------------
COPY --from=gifski-builder /usr/local/cargo/bin/gifski /usr/local/bin/gifski

# ---------------------------------------------------------------------------
# Ubuntu 22.04: GCC 13 (PPA) + newer CMake (Kitware APT)
#
# Ubuntu 24.04 ships GCC 13 and CMake 3.28+ by default — skip this step.
# ---------------------------------------------------------------------------
RUN if [ "${UBUNTU_VERSION}" = "22.04" ]; then \
      echo "==> Ubuntu 22.04 detected: installing GCC 13 PPA and Kitware CMake" \
      && add-apt-repository -y ppa:ubuntu-toolchain-r/test \
      && wget -qO- https://apt.kitware.com/keys/kitware-archive-latest.asc \
           | gpg --dearmor \
           | tee /usr/share/keyrings/kitware-archive-keyring.gpg > /dev/null \
      && echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] \
           https://apt.kitware.com/ubuntu/ jammy main" \
           | tee /etc/apt/sources.list.d/kitware.list \
      && apt-get update \
      && apt-get install -y --no-install-recommends g++-13 gcc-13 cmake \
      && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 60 \
      && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 60 \
      && rm -rf /var/lib/apt/lists/* ; \
    else \
      echo "==> Ubuntu ${UBUNTU_VERSION}: using default GCC and CMake" ; \
    fi

# ---------------------------------------------------------------------------
# Python packages for visualization scripts (auto_plot.py, etc.)
# Managed by uv from pyproject.toml + uv.lock (the single source of truth).
# The two files are copied before the full source tree so this dependency
# layer stays cached across source-only rebuilds.
# ---------------------------------------------------------------------------
COPY --from=ghcr.io/astral-sh/uv:0.10.9 /uv /usr/local/bin/uv
ENV UV_PYTHON_DOWNLOADS=never
WORKDIR /workspace
COPY pyproject.toml uv.lock ./
RUN uv sync --locked
ENV PATH="/workspace/.venv/bin:${PATH}"

# ---------------------------------------------------------------------------
# Copy project source
# ---------------------------------------------------------------------------
WORKDIR /workspace
COPY . .

# ---------------------------------------------------------------------------
# Debug build
# ---------------------------------------------------------------------------
RUN cmake \
      -G Ninja \
      -S . -B build-debug \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_COMPILER=gcc \
      -DCMAKE_CXX_COMPILER=g++ \
      -DCMAKE_CXX_STANDARD=20 \
      -DCMAKE_CXX_STANDARD_REQUIRED=ON \
      -DCMAKE_CXX_EXTENSIONS=OFF \
      -DCMAKE_C_COMPILER_LAUNCHER=ccache \
      -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
 && cmake --build build-debug -j"$(nproc)"

# ---------------------------------------------------------------------------
# Release build
# ---------------------------------------------------------------------------
RUN cmake \
      -G Ninja \
      -S . -B build-release \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=gcc \
      -DCMAKE_CXX_COMPILER=g++ \
      -DCMAKE_CXX_STANDARD=20 \
      -DCMAKE_CXX_STANDARD_REQUIRED=ON \
      -DCMAKE_CXX_EXTENSIONS=OFF \
      -DCMAKE_C_COMPILER_LAUNCHER=ccache \
      -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
      -DNODEBUG=ON \
 && cmake --build build-release -j"$(nproc)"

RUN echo 'PS1="(docker:muonith) \w# "' >> /root/.bashrc
CMD ["bash"]
