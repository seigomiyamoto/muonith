{
  description = "muonith - C++20/CMake + Python development environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };
        isLinux = pkgs.stdenv.isLinux;
        isDarwin = pkgs.stdenv.isDarwin;
      in
      {
        devShells.default = pkgs.mkShell {
          nativeBuildInputs = [
            pkgs.cmake
            pkgs.ninja
            pkgs.ccache
            pkgs.pkg-config
            pkgs.git
            pkgs.bash
          ] ++ pkgs.lib.optionals isLinux [
            pkgs.gcc
          ];

          buildInputs = pkgs.lib.optionals isLinux [
            pkgs.openblas
            pkgs.lapack
          ] ++ pkgs.lib.optionals isDarwin [
            pkgs.llvmPackages.openmp
          ];

          packages = [
            pkgs.gifski
            pkgs.imagemagick
            pkgs.doxygen
            pkgs.graphviz
            # Python packages are managed by uv (pyproject.toml + uv.lock at the
            # repo root). The flake provides only the interpreter and uv itself;
            # run "uv sync" to populate .venv.
            pkgs.python3
            pkgs.uv
          ];

          shellHook = ''
            export IN_NIX_SHELL=1
            # Use the nix-provided interpreter; never download a managed Python.
            export UV_PYTHON_DOWNLOADS=never
          '';
        };
      }
    );
}
