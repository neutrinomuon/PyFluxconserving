from numpy.distutils.core import Extension
from numpy.distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

ext1 = Extension(  name='fluxconserving.flib',
                   sources=[ 'src/fortran/DataTypes.f90',
                             'src/fortran/LINinterpol.f90',
                             'src/fortran/AkimaSpline.f90',
                             'src/fortran/SPLINE1DArr.f90',
                             'src/fortran/SPLINE3DFor.f90',
                             'src/fortran/Interpolado.f90',
                             'src/fortran/LINdexerpol.f90',
                             'src/fortran/FluxConSpec.f90' ]
                 )
    
setup( name='fluxconserving',
       version='0.0.1',
       ext_modules=[ ext1 ],
       extra_compile_args=['-O3'],
       description='Flux-conserving legacy routines in Fortran and Python',
       long_description=long_description,      # Long description read from the the readme file
       long_description_content_type="text/markdown",
       author='Jean Gomes',
       author_email='antineutrinomuon@gmail.com',
       url='https://github.com/neutrinomuon/FluxConserving',
       install_requires=[ 'numpy','matplotlib' ],
       classifiers=[
           "Programming Language :: Python :: 3",
           "Programming Language :: Fortran",
           "Operating System :: OS Independent",
                   ],
       package_dir={"fluxconserving": "src/python"},
       packages=['fluxconserving'],
       data_files=[('', ['version.txt'])],
      )
    
