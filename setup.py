'''Setup for PyFluxconserving'''
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

with open("README.md", "r",encoding="utf-8") as fh:
    long_description = fh.read()

with open("version.txt", "r",encoding="utf-8") as vh:
    version_description = vh.read()

ext1 = Extension(  name='PyFluxconserving.flib',
                   sources=[ 'src/fortran/DataTypes.f90',
                             'src/fortran/AkimaSpline.f90',
                             'src/fortran/Interpolado.f90',
                             'src/fortran/LINdexerpol.f90',
                             'src/fortran/LINinterpol.f90',
                             'src/fortran/SPLINE1DArr.f90',
                             'src/fortran/SPLINE3DFor.f90',
                             'src/fortran/FluxConSpec.f90' ]
                 )

ext2 = Extension(name='pyfluxconserving.flib',
                 sources=[ 'src/fortran/DataTypes.f90',
                           'src/fortran/AkimaSpline.f90',
                           'src/fortran/Interpolado.f90',
                           'src/fortran/LINdexerpol.f90',
                           'src/fortran/LINinterpol.f90',
                           'src/fortran/SPLINE1DArr.f90',
                           'src/fortran/SPLINE3DFor.f90',
                           'src/fortran/FluxConSpec.f90'])

setup( name='PyFluxconserving',
       version=version_description,
       ext_modules=[ ext1,ext2 ],
       extra_compile_args=['-O3'],
       description='FluxConserving is a set of Fortran 2003+ legacy routines with Python. \
There are some options for the flux-conserving algorithm. It also includes interpolation scripts.',
       long_description=long_description,      # Long description read from the the readme file
       long_description_content_type="text/markdown",
       author='Jean Gomes',
       author_email='antineutrinomuon@gmail.com',
       maintainer='Jean Gomes',
       maintainer_email='antineutrinomuon@gmail.com',
       keywords='spectra,spectral analysis,spectral synthesis,interpolation,stars,galaxies',
       url='https://github.com/neutrinomuon/PyFluxconserving',
       docs_url='https://github.com/neutrinomuon/PyFluxconserving',
       download_url='https://github.com/neutrinomuon/PyFluxconserving',
       install_requires=[ 'numpy','matplotlib','scipy' ],
       requires_python='>=3.9',
       classifiers=[
           "Programming Language :: Python :: 3",
           "Programming Language :: Fortran",
           "Operating System :: OS Independent",
                   ],
       #package_dir={"PyFluxconserving": "src/python"},
       #packages=['PyFluxconserving'],
       package_dir={"PyFluxconserving": "src/python", "pyfluxconserving": "src/python"},
       packages=['PyFluxconserving', 'pyfluxconserving'],
       data_files=[('', ['version.txt','LICENSE.txt'])],
      )
