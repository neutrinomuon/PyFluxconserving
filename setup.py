from numpy.distutils.core import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import sys
import os
import glob
import shutil
from pathlib import Path

# --- Metadata ---
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

PACKAGE_NAME = "PyFluxconserving"
MODULE_NAME = "flib"

# --- Version from version.txt ---
def read_version():
    version_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), "version.txt")
    with open(version_file, "r", encoding="utf-8") as f:
        return f.read().strip()

# --- Fortran sources ---
fortran_sources = [
    'src/fortran/DataTypes.f90',
    'src/fortran/AkimaSpline.f90',
    'src/fortran/Interpolado.f90',
    'src/fortran/LINdexerpol.f90',
    'src/fortran/LINinterpol.f90',
    'src/fortran/SPLINE1DArr.f90',
    'src/fortran/SPLINE3DFor.f90',
    'src/fortran/FluxConSpec.f90',
]

# --- Define the F2Py Extension ---
ext_modules = [
    Extension(
        name=f"{PACKAGE_NAME}.{MODULE_NAME}",
        sources=fortran_sources,
    )
]

# --- Optional custom build_ext for copying .py files ---
class F2PyBuild(build_ext):
    def run(self):
        # Ensure __init__.py exists
        pkg_dir = os.path.join(os.getcwd(), PACKAGE_NAME)
        os.makedirs(pkg_dir, exist_ok=True)
        init_file = os.path.join(pkg_dir, "__init__.py")
        if not os.path.exists(init_file):
            open(init_file, "a").close()

        # --- Generate .pyf interface using f2py ---
        cmd = [sys.executable, "-m", "numpy.f2py", "-m", MODULE_NAME, "-h",
               os.path.join(pkg_dir, f"{MODULE_NAME}.pyf")] + fortran_sources
        subprocess.check_call(cmd)

        # --- Compile Fortran to .so ---
        cmd = [sys.executable, "-m", "numpy.f2py", "-c",
               os.path.join(pkg_dir, f"{MODULE_NAME}.pyf")] + fortran_sources
        subprocess.check_call(cmd)

        super().run()

        # --- Move compiled .so to package if needed ---
        for so_file in glob.glob(f"{MODULE_NAME}*.so"):
            shutil.move(so_file, os.path.join(pkg_dir, so_file))

        # --- Copy additional Python sources into the package ---
        py_src_dir = os.path.join(os.getcwd(), "src", "python")
        if os.path.isdir(py_src_dir):
            for py_file in glob.glob(os.path.join(py_src_dir, "*.py")):
                shutil.copy(py_file, pkg_dir)

# --- Setup ---
setup(
    name=PACKAGE_NAME,
    version=read_version(),
    packages=[PACKAGE_NAME],
    cmdclass={"build_ext": F2PyBuild},
    ext_modules=ext_modules,              # Important: register extension
    has_ext_modules=lambda: True,         # Force wheel to be non-pure
    package_data={PACKAGE_NAME: ["flib*.so"]},  # Include .so files
    long_description=long_description,
    long_description_content_type="text/markdown",
)
