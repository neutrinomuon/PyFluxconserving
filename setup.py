from numpy.distutils.core import Extension
from numpy.distutils.core import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

ext1 = Extension(  name='pyfluxconserving.flib',
                   sources=[ 'src/fortran/DataTypes.f90',
                             'src/fortran/AkimaSpline.f90',
                             'src/fortran/Interpolado.f90',
                             'src/fortran/LINdexerpol.f90',
                             'src/fortran/LINinterpol.f90',
                             'src/fortran/SPLINE1DArr.f90',
                             'src/fortran/SPLINE3DFor.f90',
                             'src/fortran/FluxConSpec.f90' ]
                 )
    
setup( name='pyfluxconserving',
       version='0.0.7',
       ext_modules=[ ext1 ],
       extra_compile_args=['-O3'],
       description='FluxConserving is a set of Fortran 2003+ legacy routines with Python. There are some options for the flux-conserving algorithm. It also includes interpolation scripts.',
       long_description=long_description,      # Long description read from the the readme file
       long_description_content_type="text/markdown",
       author='Jean Gomes',
       author_email='antineutrinomuon@gmail.com',
       url='https://github.com/neutrinomuon/Pyfluxconserving',
       install_requires=[ 'numpy','matplotlib' ],
       classifiers=[
           "Programming Language :: Python :: 3",
           "Programming Language :: Fortran",
           "Operating System :: OS Independent",
                   ],
       package_dir={"pyfluxconserving": "src/python"},
       packages=['pyfluxconserving'],
       data_files=[('', ['version.txt'])],
      )
    
