from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import subprocess
import sys
import os
import glob
import shutil
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

# --- version handling from version.txt ----------------
def read_version():
    here = os.path.abspath(os.path.dirname(__file__))
    version_file = os.path.join(here, "version.txt")
    with open(version_file, "r", encoding="utf-8") as f:
        return f.read().strip()
# --- version handling from version.txt ----------------

PACKAGE_NAME = "PyFluxconserving"
MODULE_NAME = "flib"

class F2PyBuild(build_ext):
    def run(self):
        # Diretório do pacote Python
        pkg_dir = os.path.join(os.getcwd(), PACKAGE_NAME)
        os.makedirs(pkg_dir, exist_ok=True)
        init_file = os.path.join(pkg_dir, "__init__.py")
        if not os.path.exists(init_file):
            open(init_file, "a").close()  # cria __init__.py vazio

        # Lista de fontes Fortran
        sources = [
            'src/fortran/DataTypes.f90',
            'src/fortran/AkimaSpline.f90',
            'src/fortran/Interpolado.f90',
            'src/fortran/LINdexerpol.f90',
            'src/fortran/LINinterpol.f90',
            'src/fortran/SPLINE1DArr.f90',
            'src/fortran/SPLINE3DFor.f90',
            'src/fortran/FluxConSpec.f90',
        ]

        # --- Gerar .pyf interface ---
        cmd = [sys.executable, "-m", "numpy.f2py", "-m", MODULE_NAME, "-h",
               os.path.join(pkg_dir, f"{MODULE_NAME}.pyf")] + sources
        subprocess.check_call(cmd)

        # --- Compilar Fortran para .so dentro do pacote ---
        cmd = [sys.executable, "-m", "numpy.f2py", "-c",
               os.path.join(pkg_dir, f"{MODULE_NAME}.pyf")] + sources
        subprocess.check_call(cmd)

        super().run()

        # --- Mover o .so para dentro do pacote, se não estiver ---
        for so_file in glob.glob(f"{MODULE_NAME}*.so"):
            shutil.move(so_file, os.path.join(pkg_dir, so_file))


        # --- Copiar arquivos Python para dentro do pacote ---
        py_src_dir = os.path.join(os.getcwd(), "src", "python")
        if os.path.isdir(py_src_dir):
            for py_file in glob.glob(os.path.join(py_src_dir, "*.py")):
                shutil.copy(py_file, pkg_dir)

# Setup principal
setup(
    name=PACKAGE_NAME,
    version=read_version(),
    packages=[PACKAGE_NAME],
    #packages=find_packages(),  # busca pacotes em src/python
    #package_dir={"": "src/python"},              # root do pacote é src/python
    cmdclass={"build_ext": F2PyBuild},
    package_data={"PyFluxconserving": ["flib*.so"]},  # include any flib .so
    #install_requires=["numpy", "scipy", "matplotlib"],  # dependências
    long_description=long_description,
    long_description_content_type="text/markdown",
)
