# Nix Development Environment (macOS, Optional)

This guide covers setting up the MUONITH development environment using
[Nix](https://nixos.org/), a purely functional package manager that
provides **reproducible, declarative** development environments.

!!! note "Nix is optional"
    MUONITH builds natively on macOS (Apple Clang + Accelerate) and Linux
    (GCC + OpenBLAS) without Nix. See [Installation](installation.md) for
    the standard build instructions. Nix is an alternative for users who
    prefer an isolated environment that does not install packages globally.

!!! tip "Why choose Nix?"
    With Nix you get the exact same toolchain (compiler, CMake, Python
    packages, etc.) on every machine -- no manual dependency installation,
    no version mismatches, and no "works on my machine" surprises.
    A single `nix develop` command sets up everything without modifying
    your system.

## Prerequisites

### Install Nix

The recommended installer is the
[Determinate Systems installer](https://github.com/DeterminateSystems/nix-installer),
which enables flakes by default and provides a clean uninstall path.

```bash
curl --proto '=https' --tlsv1.2 -sSf -L \
  https://install.determinate.systems/nix | sh -s -- install
```

After installation, open a new terminal (or `source` your shell profile)
and verify:

```bash
nix --version
```

!!! warning "Flakes must be enabled"
    The Determinate Systems installer enables flakes automatically.
    If you installed Nix via the official installer instead, you need to
    enable flakes manually. Add the following to `~/.config/nix/nix.conf`:

    ```ini
    experimental-features = nix-command flakes
    ```

## Getting started

### Enter the development shell

```bash
cd ~/prg/muonith       # Replace with the actual path to your clone
nix develop
```

The first invocation downloads and builds all dependencies (this may
take several minutes). Subsequent invocations use the Nix store cache
and are near-instant.

### Verify the environment

Once inside the shell, confirm that the key tools are available:

```bash
cmake --version
ninja --version
uv --version
echo $IN_NIX_SHELL   # Should print "1"
```

Python packages are not part of the flake; set them up once with:

```bash
uv sync                  # creates .venv from pyproject.toml + uv.lock
.venv/bin/python -c "import numpy; print(numpy.__version__)"
```

## Platform differences

The `flake.nix` adapts to the host platform automatically. The key
differences are:

| | macOS | Linux |
|---|---|---|
| **Compiler** | Clang (via Nix stdenv) | GCC |
| **BLAS / LAPACK** | Apple Accelerate | OpenBLAS + LAPACK |
| **OpenMP** | `llvmPackages.openmp` | Included with GCC |
| **Nix wrappers** | `cc` = Clang, `c++` = Clang++ | `cc` = GCC, `c++` = G++ |

!!! note "No Homebrew GCC required"
    Unlike the [manual installation](installation.md) path, the Nix
    environment does **not** require Homebrew GCC on macOS. Clang from
    Nix stdenv provides full C++20 support and OpenMP is supplied as a
    separate `llvmPackages.openmp` package.

## Building MUONITH

Inside the Nix shell, use the standard build scripts:

```bash
# Debug build
bash build.sh Debug

# Or use the wrapper script
bash bdebug.sh
```

```bash
# Release build
bash build.sh Release

# Or use the wrapper script
bash brelease.sh
```

The build scripts detect the `IN_NIX_SHELL` environment variable and
select the appropriate compiler automatically -- no manual
`CMAKE_C_COMPILER` or `CMAKE_CXX_COMPILER` flags needed.

## The `IN_NIX_SHELL` environment variable

The `flake.nix` shellHook exports `IN_NIX_SHELL=1` when you enter the
development shell. The `build.sh` script uses this variable to decide
compiler selection:

| Condition | `CC` | `CXX` |
|---|---|---|
| `IN_NIX_SHELL` is set | `cc` (Nix wrapper) | `c++` (Nix wrapper) |
| macOS without Nix | `gcc-14` or `gcc-13` (Homebrew) | `g++-14` or `g++-13` |
| Linux without Nix | `gcc` | `g++` |

This means the same `build.sh` works identically inside and outside of
Nix, with no manual configuration.

## What the flake provides

For reference, the complete set of packages available in the Nix
development shell:

| Category | Packages |
|---|---|
| **Build tools** | CMake, Ninja, ccache, pkg-config, git, bash |
| **Compiler** | GCC (Linux only; macOS uses Clang via stdenv) |
| **Libraries** | OpenBLAS + LAPACK (Linux), llvmPackages.openmp (macOS) |
| **Documentation** | Doxygen, Graphviz |
| **Media** | gifski (GIF encoding), ImageMagick (`convert`, for PDF/GIF figure output) |
| **Python** | Python 3 interpreter and `uv` (no packages) |

Python packages are intentionally not provided by the flake. They are
managed by `uv sync` from `pyproject.toml` + `uv.lock` at the repository
root, so every platform (with or without Nix) installs the same pinned
set. The shellHook sets `UV_PYTHON_DOWNLOADS=never` so that `uv` always
uses the Nix-provided interpreter instead of downloading its own.

## Troubleshooting

### Flakes not enabled

If you see an error like:

```
error: experimental Nix feature 'flakes' is disabled
```

Add the following to `~/.config/nix/nix.conf` (create the file if it
does not exist):

```ini
experimental-features = nix-command flakes
```

Then restart your terminal.

!!! tip
    If you used the [Determinate Systems installer](https://github.com/DeterminateSystems/nix-installer),
    flakes are enabled by default and you should not see this error.

### Slow first `nix develop`

The first invocation of `nix develop` downloads and builds all
dependencies from source or binary cache. This can take several minutes
depending on your network speed.

Subsequent invocations are fast because Nix caches everything in the
local store (`/nix/store/`).

!!! tip "Binary cache"
    Nix automatically uses the official binary cache
    (`cache.nixos.org`). Most packages are pre-built and download as
    binaries. If a specific package is not cached (e.g., a pinned Python
    package), Nix builds it locally.

### `ccache` not finding the cache

Inside the Nix shell, `ccache` uses its default cache directory
(`~/.cache/ccache` or `~/.ccache`). If you have a non-Nix ccache
installation with a different cache directory, the caches are shared
automatically.

Verify ccache status:

```bash
ccache -s
```

### macOS: OpenMP not found

If CMake reports that OpenMP cannot be found on macOS, ensure you are
inside the Nix shell (`echo $IN_NIX_SHELL` should print `1`). The Nix
shell provides `llvmPackages.openmp` and sets the necessary paths
automatically.

If building **outside** Nix on macOS, see the
[manual installation guide](installation.md#openmp) instead.
