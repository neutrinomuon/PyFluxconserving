# Example of compilation
# LAPACK library is in use

f2py -c --opt=-O3 --opt=-lapack ../Modules/DataTypes.f90 AkimaSpline.f90 LINinterpol.f90 pchipmodule.f90 SPLINE1DArr.f90 SPLINE3DFor.f90 Interpolado.f90 LINdexerpol.f90 Fitting_Fun.f90 FluxConSpec.f90 -m FluxConSpec

