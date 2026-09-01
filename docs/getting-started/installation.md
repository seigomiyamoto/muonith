# Installation

## Quick start — just run these commands

Copy-paste the commands for your platform to go from zero to a working build.

=== "Ubuntu 24.04"

    === "native"

        ```bash
        # 1. System packages
        sudo apt update
        sudo apt install -y \
          build-essential cmake ninja-build ccache pkg-config git \
          libopenblas-openmp-dev liblapacke-dev \
          time util-linux coreutils python3-venv \
          imagemagick

        # 2. Clone and enter the repository
        git clone git@github.com:seigomiyamoto/muonith.git   # SSH
        # git clone https://github.com/seigomiyamoto/muonith.git  # HTTPS
        cd muonith

        # 3. Build
        bash bdebug.sh      # Debug build
        bash brelease.sh    # Release build (optimized)

        # 4. Python packages for visualization (optional)
        curl -LsSf https://astral.sh/uv/install.sh | sh
        source ~/.local/bin/env  # or restart your shell
        uv sync                  # creates .venv from pyproject.toml + uv.lock
        source .venv/bin/activate
        ```

    === "Docker"

        ```bash
        # 1. Install Docker
        sudo apt update
        sudo apt install -y docker.io git
        sudo usermod -aG docker $USER
        ```

        Log out and back in, then:

        ```bash
        # 2. Clone and build the Docker image
        git clone git@github.com:seigomiyamoto/muonith.git   # SSH
        # git clone https://github.com/seigomiyamoto/muonith.git  # HTTPS
        cd muonith
        docker build -f docker/Dockerfile -t muonith-test .

        # 3. Verify
        docker run -it muonith-test bash
        # Inside: ls build-debug/exec/*.exe && ls build-release/exec/*.exe && exit
        ```

        See [Docker details](#alternative-build-with-docker-linux-wsl2) for
        VS Code integration and cleanup.

=== "WSL2 (Ubuntu 24.04)"

    === "native"

        ```bash
        # 1. System packages
        sudo apt update
        sudo apt install -y \
          build-essential cmake ninja-build ccache pkg-config git \
          libopenblas-openmp-dev liblapacke-dev \
          time util-linux coreutils python3-venv \
          imagemagick

        # 2. Clone and enter the repository (use /home/, NOT /mnt/c/)
        cd ~
        git clone git@github.com:seigomiyamoto/muonith.git   # SSH
        # git clone https://github.com/seigomiyamoto/muonith.git  # HTTPS
        cd muonith

        # 3. Build
        bash bdebug.sh      # Debug build
        bash brelease.sh    # Release build (optimized)

        # 4. Python packages for visualization (optional)
        curl -LsSf https://astral.sh/uv/install.sh | sh
        source ~/.local/bin/env  # or restart your shell
        uv sync                  # creates .venv from pyproject.toml + uv.lock
        source .venv/bin/activate
        ```

        !!! warning "Store files under `/home/`, not `/mnt/c/`"
            The Windows filesystem (`/mnt/c/`) is 10–15× slower. Always work
            under the native Linux filesystem.

    === "Docker"

        ```bash
        # 1. Install Docker
        sudo apt update
        sudo apt install -y docker.io git
        sudo usermod -aG docker $USER
        ```

        Restart WSL2 (`wsl --shutdown` in PowerShell), then reopen the terminal:

        ```bash
        # 2. Clone and build the Docker image (use /home/, NOT /mnt/c/)
        cd ~
        git clone git@github.com:seigomiyamoto/muonith.git   # SSH
        # git clone https://github.com/seigomiyamoto/muonith.git  # HTTPS
        cd muonith
        docker build -f docker/Dockerfile -t muonith-test .

        # 3. Verify
        docker run -it muonith-test bash
        # Inside: ls build-debug/exec/*.exe && ls build-release/exec/*.exe && exit
        ```

        See [Docker details](#alternative-build-with-docker-linux-wsl2) for
        VS Code integration and cleanup.

=== "macOS"

    === "Nix"

        Install Nix first if you haven't — see [Nix setup guide](nix-setup.md#install-nix).

        ```bash
        # 1. Clone and enter the repository
        git clone git@github.com:seigomiyamoto/muonith.git   # SSH
        # git clone https://github.com/seigomiyamoto/muonith.git  # HTTPS
        cd muonith

        # 2. Enter Nix dev shell (provides compiler, cmake, etc.)
        nix develop

        # 3. Build
        bash bdebug.sh      # Debug build
        bash brelease.sh    # Release build (optimized)

        # 4. Python packages for visualization (optional)
        uv sync                  # creates .venv from pyproject.toml + uv.lock
        source .venv/bin/activate
        ```

    === "native"

        ```bash
        # 1. Xcode CLI Tools + system packages
        xcode-select --install
        brew install cmake ninja ccache pkg-config git gnu-time libomp gifski imagemagick

        # 2. Clone and enter the repository
        git clone git@github.com:seigomiyamoto/muonith.git   # SSH
        # git clone https://github.com/seigomiyamoto/muonith.git  # HTTPS
        cd muonith

        # 3. Build
        bash bdebug.sh      # Debug build
        bash brelease.sh    # Release build (optimized)

        # 4. Python packages for visualization (optional)
        curl -LsSf https://astral.sh/uv/install.sh | sh
        source ~/.local/bin/env  # or restart your shell
        # Or: brew install uv
        uv sync                  # creates .venv from pyproject.toml + uv.lock
        source .venv/bin/activate
        ```

For detailed explanations of each dependency, read on below.

---

This guide covers building MUONITH from source on the following platforms:

| Platform | Status |
|---|---|
| Linux (Ubuntu 22.04 / 24.04, including WSL2) | Fully supported |
| Linux (Docker) | Fully supported |
| macOS (Apple Clang + Accelerate) | Fully supported |
| macOS (Nix development environment) | Fully supported |

## Choose your build method

| Platform | Recommended | Alternative |
|---|---|---|
| **Linux / WSL2** | [Docker](#alternative-build-with-docker-linux-wsl2) (reproducible, no host modification) | Install packages via `apt` (this page) |
| **macOS** | [Nix](nix-setup.md) (isolated, reproducible environment) | Install via Homebrew + Xcode CLI Tools (this page) |

!!! success "macOS native build verified"
    macOS native builds with Apple Clang 17 + Accelerate framework have been
    verified on macOS 15 Sequoia (Apple Silicon M4 Pro). All 110 build targets
    pass and runtime tests complete successfully. No Homebrew GCC or OpenBLAS
    is required.

The sections below describe the **native** (alternative) installation:

- **Linux / WSL2 users** — [Docker](#alternative-build-with-docker-linux-wsl2) is recommended
  for reproducibility and zero host modification. The native `apt` instructions
  below are provided for users who prefer direct installation.
- **macOS users** — [Nix development environment](nix-setup.md) is recommended.
  It provides all dependencies in an isolated shell without installing
  Homebrew packages globally. The native instructions below are an alternative.

---

## Prerequisites

| Dependency | Purpose | Install |
|---|---|---|
| **C++20 compiler** | GCC 13+ or Apple Clang 17+ | See [Compiler](#compiler) below |
| **CMake** | Build system | 3.22.2+ |
| **Ninja** | Build generator | — |
| **BLAS/LAPACK** | Linear algebra | Platform-dependent (see [below](#blaslapack)) |
| **OpenMP** | Thread parallelization | (included with GCC; separate on macOS) |
| **ccache** *(optional)* | Compiler cache for faster rebuilds | — |
| **Doxygen** *(optional)* | API documentation generation | — |

!!! note "C++ libraries downloaded automatically"
    **Eigen3**, **fmt**, **spdlog**, and **nlohmann/json** are automatically
    downloaded and built by CMake via FetchContent. No manual installation is
    required.

## System packages

=== "Linux / WSL2"

    ```bash
    sudo apt update
    sudo apt install -y \
      build-essential \
      cmake \
      ninja-build \
      ccache \
      pkg-config \
      git \
      time \
      python3-venv \
      imagemagick
    ```

    `time` provides `/usr/bin/time -v` (used by station `run_prg.sh`
    scripts; the bash builtin `time` is not equivalent). `python3-venv`
    is required by `setup_station.sh` to auto-create a `.venv/` for the
    visualization stack. `imagemagick` is used by `auto_plot.py` for
    PDF figure generation. Animated-GIF output additionally needs
    `gifski`, which is not packaged in `apt` — install it separately from
    its official prebuilt static binary
    (see [Visualization binaries](#visualization-binaries-optional) below).

    ??? tip "WSL2-specific tips"
        **Performance**: Always store your project files under `/home/` (the
        native Linux filesystem), not `/mnt/c/` (the Windows filesystem).
        The Windows filesystem is 10–15× slower for build operations.

        **Memory**: By default WSL2 uses 50 % of host RAM. If you run into
        memory issues during the build, create or edit
        `%USERPROFILE%\.wslconfig` on the Windows side:

        ```ini
        [wsl2]
        memory=16GB
        ```

        Then restart WSL2 with `wsl --shutdown`.

=== "macOS"

    Install [Homebrew](https://brew.sh/) if not already present, then:

    ```bash
    brew install cmake ninja ccache pkg-config git gnu-time gifski imagemagick
    ```

    `gnu-time` provides the `gtime` command that `run_prg.sh` uses on macOS
    to run the solver with time/memory measurement.

    `gifski` and `imagemagick` are optional — needed only for animated-GIF
    and PDF figure generation in `auto_plot.py`. See
    [Visualization binaries](#visualization-binaries-optional) below.

## Compiler

MUONITH requires a C++20-capable compiler. The build script `build.sh` automatically
selects the compiler based on the environment:

| Environment | Compiler | Selection logic |
|---|---|---|
| macOS (native) | Apple Clang (`/usr/bin/clang++`) | Default; override with `CC`/`CXX` |
| macOS (Nix shell) | `cc` / `c++` (Nix wrapper) | Apple Clang via stdenv |
| Linux (native) | `gcc` / `g++` | System default |
| Linux (Nix shell) | `cc` / `c++` (Nix wrapper) | GCC via stdenv |

!!! note "Apple Clang on macOS"
    Apple Clang 17+ (Xcode Command Line Tools) fully supports the C++20 features
    used by MUONITH and is what `build.sh` selects on native macOS. The absolute
    path is used so a Homebrew LLVM clang earlier on `PATH` is not picked up.
    OpenMP requires `libomp` (see [OpenMP](#openmp)).

=== "Linux / WSL2 (Ubuntu 24.04)"

    G++ 13 is included by default with Ubuntu 24.04. No additional steps needed.

    Verify:

    ```bash
    g++ --version
    # g++ (Ubuntu 13.3.0-...) 13.3.0
    ```

=== "Linux / WSL2 (Ubuntu 22.04)"

    ```bash
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test
    sudo apt update
    sudo apt install g++-13 gcc-13
    ```

    Optionally set as default:

    ```bash
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 60
    sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 60
    ```

    Verify:

    ```bash
    g++ --version
    # g++ (Ubuntu 13.x.x) 13.x.x
    ```

=== "macOS"

    **Apple Clang** (included with Xcode Command Line Tools) is the default
    compiler on macOS. Install the tools if not already present:

    ```bash
    xcode-select --install
    ```

    Verify:

    ```bash
    clang++ --version
    # Apple clang version 17.0.0 (...)
    ```

    Apple Clang 17+ fully supports the C++20 features used by MUONITH.
    OpenMP is **not** bundled with Apple Clang — install `libomp` separately
    (see [OpenMP](#openmp)).

    !!! note "Optional: Homebrew GCC"
        Homebrew GCC is not selected automatically -- the native macOS build is
        verified only with Apple Clang. To use it anyway, install it with
        `brew install gcc@14` and export the compilers explicitly:
        `CC=gcc-14 CXX=g++-14 bash bdebug.sh`. GCC bundles OpenMP, so no
        separate `libomp` is needed in that case.

## BLAS/LAPACK

MUONITH uses BLAS/LAPACK for matrix operations in the inversion step. The backend
is selected automatically by CMake based on the platform:

| Platform | BLAS/LAPACK backend | Notes |
|---|---|---|
| macOS | Apple Accelerate (default) | Built into macOS, no installation needed |
| macOS (optional) | Homebrew OpenBLAS | Requires `CMAKE_PREFIX_PATH` |
| Linux | OpenBLAS + LAPACK | Manual installation required |

=== "Linux / WSL2"

    ```bash
    sudo apt install -y libopenblas-openmp-dev liblapacke-dev
    ```

    Verify:

    ```bash
    dpkg -l | grep libopenblas
    ```

=== "macOS"

    **Apple Accelerate** is built into macOS — no installation is required.
    CMake automatically detects and links the Accelerate framework.

    !!! note "Nix users"
        In the Nix shell, macOS also uses Apple Accelerate (via the system
        framework). No additional Nix packages are needed for BLAS/LAPACK.

!!! info "Numerical differences between backends"
    Accelerate and OpenBLAS use different internal algorithms. This can produce
    tiny numerical differences (symmetric relative error ~10⁻⁷) in density
    reconstruction results. Both are correct — see [Inversion: Numerical
    Implementation](../concepts/inversion.md#numerical-implementation) for details.

## OpenMP

=== "Linux / WSL2"

    OpenMP is included with GCC — no additional installation is required.

=== "macOS"

    If compiling with **Apple Clang** (the default), install `libomp`:

    ```bash
    brew install libomp
    ```

    `build.sh` automatically detects Homebrew libomp and passes
    `-DOpenMP_ROOT` to CMake.

    If you override the compiler with a Homebrew GCC (`CC=gcc-14`), OpenMP is
    included with it — no additional package is needed.

## Doxygen (optional)

Only needed if you want to generate API reference documentation.

=== "Linux / WSL2"

    ```bash
    sudo apt install doxygen graphviz
    ```

=== "macOS"

    ```bash
    brew install doxygen graphviz
    ```

## Python packages (optional)

Python 3 is required for visualization and post-processing scripts
(`auto_plot.py`, `hist2d.py`, `plot_g2bg.py`, etc.).
These are **not needed** to build or run MUONITH itself.

!!! tip "Recommended: use `uv` for virtual environment management"
    [`uv`](https://docs.astral.sh/uv/) is a fast Python package manager
    written in Rust. It creates isolated virtual environments without
    polluting your system Python.

### Install uv

=== "Linux / WSL2"

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

=== "macOS"

    ```bash
    curl -LsSf https://astral.sh/uv/install.sh | sh
    ```

    Or via Homebrew:

    ```bash
    brew install uv
    ```

### Create a virtual environment and install packages

The package list is maintained in `pyproject.toml` at the repository root, and
exact versions are pinned in `uv.lock`. A single command reproduces the
environment:

```bash
uv sync                  # creates .venv from pyproject.toml + uv.lock
source .venv/bin/activate
```

!!! note "macOS"
    Pre-built binary packages for `rasterio` and `pyproj` require macOS 14
    (Sonoma) or newer. On older macOS versions `uv sync` falls back to
    building them from source, which needs additional system libraries.

??? note "Alternative: traditional pip (not recommended)"
    If you prefer not to install `uv`, you can use the system `pip`
    directly. Note that installing packages globally with `pip` is
    [deprecated in recent Python versions](https://peps.python.org/pep-0668/)
    and may require `--break-system-packages` on Ubuntu 24.04+.

    === "Linux / WSL2"

        ```bash
        sudo apt install -y python3 python3-pip python3-venv
        python3 -m venv .venv
        source .venv/bin/activate
        pip install numpy matplotlib pandas scipy requests json5 pyproj rasterio pillow
        ```

    === "macOS"

        ```bash
        python3 -m venv .venv
        source .venv/bin/activate
        pip install numpy matplotlib pandas scipy requests json5 pyproj rasterio pillow
        ```

### Documentation development (mkdocs)

To build or preview this documentation site locally:

```bash
uv sync --group docs     # adds the mkdocs toolchain to .venv
source .venv/bin/activate
mkdocs serve             # Preview at http://127.0.0.1:8001
```

!!! tip "Additional packages for GIS tools"
    The auxiliary tools `muonith-path-view` and `muonith-gsi-dem` require extra
    packages. See [muonith-path-view](../auxiliary-tools/path-view.md) and
    [muonith-gsi-dem](../auxiliary-tools/gsi-dem.md) for details.

## Visualization binaries (optional)

The post-processing visualization driver `auto_plot.py` (invoked at the end
of `run_prg.sh`) shells out to two external binaries:

| Binary | Purpose | Used by |
|---|---|---|
| `gifski` | Combine PNG frames into animated GIFs of cross-section sweeps | `scripts/plot_Grid3dVoxel.py:638-740` (`run_gifski_cross`) |
| `convert` (ImageMagick) | Combine PNG frames into multi-page PDFs of cross-section sweeps | `scripts/plot_Grid3dVoxel.py:727-744` (`run_convert_pdf_cross`) |

The PNG figures themselves are always generated by matplotlib and do **not**
require these binaries. Only the animated GIF and PDF outputs require them.

=== "Linux / WSL2"

    `imagemagick` (the `convert` command) is available via `apt`, but
    `gifski` is **not** packaged in Ubuntu's `apt`. Install `imagemagick`
    from `apt`, then install `gifski` from its official prebuilt static
    Linux binary — no Rust toolchain required:

    ```bash
    sudo apt install -y imagemagick curl xz-utils
    mkdir -p ~/.local/bin
    curl -sL -o /tmp/gifski.tar.xz \
      https://github.com/ImageOptim/gifski/releases/download/1.34.0/gifski-1.34.0.tar.xz
    tar -xf /tmp/gifski.tar.xz -C /tmp linux/gifski
    install -m 755 /tmp/linux/gifski ~/.local/bin/gifski
    ```

    The archive bundles a fully static binary (`linux/gifski`) with no
    shared-library dependencies; it runs unchanged on Ubuntu 20.04, 22.04,
    and 24.04. The download URL is version-pinned — check the
    [gifski releases page](https://github.com/ImageOptim/gifski/releases)
    for a newer version and adjust the URL if needed. Ensure `~/.local/bin`
    is on your `PATH` (open a new shell, or
    `export PATH="$HOME/.local/bin:$PATH"`).

    Verify:

    ```bash
    which gifski && gifski --version
    which convert && convert --version | head -1
    ```

=== "macOS"

    ```bash
    brew install gifski imagemagick
    ```

    Verify:

    ```bash
    which gifski && gifski --version
    which magick && magick --version | head -1
    ```

    !!! note "ImageMagick 7 on macOS"
        Homebrew installs ImageMagick 7, which provides the `magick` command.
        A compatibility symlink for `convert` is also created by Homebrew, so
        `plot_Grid3dVoxel.py`'s direct `convert` invocation still works.

??? warning "Without these binaries"
    `auto_plot.py` will still produce all PNG figures, but the final
    GIF/PDF combination step will fail with:

    ```
    FileNotFoundError: [Errno 2] No such file or directory: 'gifski'
    ```

    (and similarly for `convert`). To recover after installing the
    binaries, re-run only the visualization step from any station
    directory:

    ```bash
    cd work/<station_name>/swp001
    bash run_prg.sh false 42   # RUN_EXE=false skips the C++ compute step
    ```

??? note "Ubuntu ImageMagick PDF policy"
    On some Ubuntu installations, ImageMagick's default
    `/etc/ImageMagick-6/policy.xml` blocks PDF output with:

    ```xml
    <policy domain="coder" rights="none" pattern="PDF" />
    ```

    If PDF generation fails with a policy error, edit that line to
    `rights="read|write"` (requires `sudo`).

## Building MUONITH

After installing all dependencies:

=== "Linux / WSL2"

    ```bash
    # Debug build
    bash bdebug.sh

    # Release build
    bash brelease.sh

    # Clean debug build (removes build directory first)
    bash ccdebug.sh
    ```

=== "macOS"

    The build scripts fully support macOS. `build.sh` detects CPU count
    via `sysctl` and selects the appropriate compiler automatically.

    ```bash
    # Debug build
    bash bdebug.sh

    # Release build
    bash brelease.sh
    ```

    !!! info "Verified on"
        macOS 15 Sequoia, Apple Silicon (M4 Pro), Apple Clang 17 + Accelerate.
        All 110 build targets pass. See also the
        [Nix development environment](nix-setup.md) as an alternative
        isolated setup.

## Alternative: Build with Docker (Linux / WSL2)

Instead of installing all dependencies on your host, you can use Docker to
build MUONITH in an isolated Ubuntu container. This is useful for:

- **Quick verification** — test that the build works without modifying your system
- **Reproducibility** — identical environment regardless of host configuration

!!! warning "macOS Docker builds not supported"
    Docker on macOS (Apple Silicon) requires Rosetta 2 emulation for the
    `linux/amd64` image, and ARM-native builds have unresolved OpenBLAS
    version compatibility issues. macOS users should use the
    [native build](#building-muonith) or
    [Nix development environment](nix-setup.md) instead.

### Install Docker

!!! note "Docker address pools may conflict with your LAN or VPN"
    In some environments, Docker's default private address pools
    (for example `172.17.0.0/16`) may overlap with your LAN or VPN routes
    and cause connectivity issues.

    If you encounter routing failures after starting Docker, reconfigure
    Docker to use a non-overlapping subnet:

    - **Linux-native Docker** — edit `/etc/docker/daemon.json` and set
      `bip` (default bridge) or `default-address-pools` (all auto-created
      networks).
    - **Docker Desktop (Windows / WSL2)** — open
      *Settings → Docker Engine* and add the same keys to the JSON config.

    Example (`daemon.json`):

    ```json
    {
      "default-address-pools": [
        {"base": "192.168.200.0/24", "size": 24}
      ]
    }
    ```

    Choose an address range that does not overlap with any subnet on your
    LAN or VPN.

=== "Linux-native"

    ```bash
    sudo apt update
    sudo apt install -y docker.io
    sudo usermod -aG docker $USER
    ```

    Log out and back in for the group change to take effect.

=== "WSL2"

    ```bash
    sudo apt update
    sudo apt install -y docker.io
    sudo usermod -aG docker $USER
    ```

    Restart WSL2 for the group change to take effect:

    ```powershell
    # In Windows PowerShell
    wsl --shutdown
    ```

    Then reopen your WSL2 terminal.

### Build the image

```bash
cd ~/prg/muonith
docker build -f docker/Dockerfile -t muonith-test .

# Ubuntu 22.04 variant (default is 24.04):
# docker build -f docker/Dockerfile --build-arg UBUNTU_VERSION=22.04 -t muonith-test .
```

First build takes 20–30 minutes (FetchContent downloads + compilation).
Subsequent builds use Docker layer cache and are much faster.

### Verify the build

```bash
# Enter the container interactively
docker run -it muonith-test bash
```

Inside the container:

```
/workspace# ls build-debug/exec/*.exe
/workspace# ls build-release/exec/*.exe
/workspace# exit
```

### Browse files with VS Code

You can attach VS Code to the running container for a full IDE experience:

1. Install the **Dev Containers** extension (`ms-vscode-remote.remote-containers`)
2. Start the container in the background:
   ```bash
   docker run -it -d --name muonith muonith-test bash
   ```
3. In VS Code: `Ctrl+Shift+P` → **Dev Containers: Attach to Running Container** → select `muonith`
4. Open folder `/workspace` — you can browse source code and build output

When finished:

```bash
docker stop muonith
docker rm muonith
```

### Clean up

```bash
# Remove the image
docker rmi muonith-test

# Remove stopped containers
docker container prune

# Remove unused images
docker image prune
```

## Generating API documentation

```bash
bash rebuild-api-docs.sh
```

HTML output is generated in `apidocs/html/`.
