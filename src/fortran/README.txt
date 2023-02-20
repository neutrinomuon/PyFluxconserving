# Example of compilation
# LAPACK library is not anymore in use since is dependend on previous installations

f2py -c --opt=-O3 DataTypes.f90 AkimaSpline.f90 LINinterpol.f90 Interpolado.f90 LINdexerpol.f90 SPLINE1DArr.f90 SPLINE3DFor.f90 FluxConSpec.f90 -m FluxConSpec


