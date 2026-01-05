name: Build Ubuntu 22.04 Package (Python 3.11 & 3.12)

on:
  push:
    branches: [ main ]
  pull_request:
    branches: [ main ]

jobs:
  build-ubuntu:
    runs-on: ubuntu-22.04
    env:
      FFLAGS: "-O2 -fallow-argument-mismatch"

    steps:
      - uses: actions/checkout@v4

      - name: Install system dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y gfortran build-essential ninja-build
          ninja --version

      - name: Build for Python 3.11
        uses: actions/setup-python@v4
        with:
          python-version: 3.11
      - run: |
          python -m pip install --upgrade pip setuptools wheel numpy meson
          # --- Compila os arquivos Fortran com f2py ---
          mkdir -p build/ubuntu-22.04-py3.11
          f2py -c src/fortran/DataTypes.f90 \
                src/fortran/AkimaSpline.f90 \
                src/fortran/Interpolado.f90 \
                src/fortran/LINdexerpol.f90 \
                src/fortran/LINinterpol.f90 \
                src/fortran/SPLINE1DArr.f90 \
                src/fortran/SPLINE3DFor.f90 \
                src/fortran/FluxConSpec.f90 \
                -m flib -h PyFluxconserving/flib.pyf
          f2py -c PyFluxconserving/flib.pyf \
                src/fortran/DataTypes.f90 \
                src/fortran/AkimaSpline.f90 \
                src/fortran/Interpolado.f90 \
                src/fortran/LINdexerpol.f90 \
                src/fortran/LINinterpol.f90 \
                src/fortran/SPLINE1DArr.f90 \
                src/fortran/SPLINE3DFor.f90 \
                src/fortran/FluxConSpec.f90 \
                -m PyFluxconserving/flib
          # --- Copia arquivos Python ---
          mkdir -p PyFluxconserving
          cp src/python/*.py PyFluxconserving/
          # --- Gera sdist + wheel no diretório correto ---
          python setup.py sdist --dist-dir dist/ubuntu-22.04-py3.11
          python setup.py bdist_wheel --dist-dir dist/ubuntu-22.04-py3.11
          # --- Move build artifacts ---
          mkdir -p build/ubuntu-22.04-py3.11
          mv build/* build/ubuntu-22.04-py3.11/ || true

      - name: Build for Python 3.12
        uses: actions/setup-python@v4
        with:
          python-version: 3.12
      - run: |
          python -m pip install --upgrade pip setuptools wheel numpy meson
          # --- Compila os arquivos Fortran com f2py ---
          mkdir -p build/ubuntu-22.04-py3.12
          f2py -c src/fortran/DataTypes.f90 \
                src/fortran/AkimaSpline.f90 \
                src/fortran/Interpolado.f90 \
                src/fortran/LINdexerpol.f90 \
                src/fortran/LINinterpol.f90 \
                src/fortran/SPLINE1DArr.f90 \
                src/fortran/SPLINE3DFor.f90 \
                src/fortran/FluxConSpec.f90 \
                -m flib -h PyFluxconserving/flib.pyf
          f2py -c PyFluxconserving/flib.pyf \
                src/fortran/DataTypes.f90 \
                src/fortran/AkimaSpline.f90 \
                src/fortran/Interpolado.f90 \
                src/fortran/LINdexerpol.f90 \
                src/fortran/LINinterpol.f90 \
                src/fortran/SPLINE1DArr.f90 \
                src/fortran/SPLINE3DFor.f90 \
                src/fortran/FluxConSpec.f90 \
                -m PyFluxconserving/flib
          # --- Copia arquivos Python ---
          mkdir -p PyFluxconserving
          cp src/python/*.py PyFluxconserving/
          # --- Gera sdist + wheel no diretório correto ---
          python setup.py sdist --dist-dir dist/ubuntu-22.04-py3.12
          python setup.py bdist_wheel --dist-dir dist/ubuntu-22.04-py3.12
          # --- Move build artifacts ---
          mkdir -p build/ubuntu-22.04-py3.12
          mv build/* build/ubuntu-22.04-py3.12/ || true

      - name: Commit & push generated files
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "github-actions[bot]@users.noreply.github.com"
          git add build dist
          git commit -m "Add build and dist for Ubuntu 22.04 Python 3.11 & 3.12" || echo "No changes to commit"
          git pull --rebase origin main || true
          git push origin main || true

      - name: Upload artifacts
        uses: actions/upload-artifact@v4
        with:
          name: all-dist
          path: |
            build/ubuntu-22.04-*
            dist/ubuntu-22.04-*
