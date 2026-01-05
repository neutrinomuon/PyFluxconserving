from setuptools import setup
import os
import glob
import shutil

PACKAGE_NAME = "PyFluxconserving"
MODULE_NAME = "flib"

# --- Lê versão do arquivo version.txt ---
def read_version():
    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "version.txt")
    with open(version_file, "r", encoding="utf-8") as f:
        return f.read().strip()

# --- Garante que o pacote exista ---
pkg_dir = os.path.join(os.getcwd(), PACKAGE_NAME)
os.makedirs(pkg_dir, exist_ok=True)

# --- Copia os .so compilados pelo CI para dentro do pacote ---
for so_file in glob.glob(f"build/{MODULE_NAME}*.so"):
    shutil.copy(so_file, pkg_dir)

# --- Copia arquivos Python do src/python para dentro do pacote ---
py_src_dir = os.path.join(os.getcwd(), "src/python")
if os.path.isdir(py_src_dir):
    for py_file in glob.glob(os.path.join(py_src_dir, "*.py")):
        shutil.copy(py_file, pkg_dir)

# --- Garante que __init__.py exista ---
init_file = os.path.join(pkg_dir, "__init__.py")
if not os.path.exists(init_file):
    open(init_file, "a").close()

# --- Setup principal ---
setup(
    name=PACKAGE_NAME,
    version=read_version(),
    packages=[PACKAGE_NAME],
    include_package_data=True,
    package_data={PACKAGE_NAME: ["flib*.so"]},  # inclui os .so compilados
    has_ext_modules=lambda: True,               # marca como pacote não-puro
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
)
