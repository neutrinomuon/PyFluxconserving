# FluxConserving
##  A fortran legacy package to easy compute the flux-density conservation
email: antineutrinomuon@gmail.com, jean@astro.up.pt

© Copyright ®

J.G. - Jean Gomes

<hr>

[![My Skills](https://skillicons.dev/icons?i=python,fortran,c,numpy&theme=light)](https://skillicons.dev)<br>
[![python3](https://img.shields.io/pypi/pyversions/pyfluxconserving)](https://img.shields.io/pypi/pyversions/pyfluxconserving)
[![badgetlicense](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/license.svg)](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/license.svg)

<hr>

<div align="center">
<img src='tutorials/fluxconserving.png' width="60%">
</div>

<br>
RESUME : Original Fortran 2003+ routines date back to 2003-2004. Read the <a
href='https://github.com/neutrinomuon/FluxConserving/blob/main/LICENSE.txt'>LICENSE.txt</a>
file. When analyzing astronomical spectra, astronomers often bin the data to
increase the signal-to-noise ratio and reduce the effects of noise in the
data. Binning refers to the process of averaging the intensity of adjacent
spectral channels, or pixels, to produce a new, coarser set of data.

In the process of binning, it is important to ensure that the principle of
flux density conservation is maintained. This means that the total energy
emitted by the object, and hence its flux density, must remain constant after
binning.

To conserve flux density, the intensity of each binned pixel should be scaled
by the number of pixels it represents. For example, if two adjacent pixels are
binned together, the intensity of the resulting bin should be the sum of the
intensities of the two original pixels, divided by two. This ensures that the
total energy in the bin is conserved, and that the flux density of the object
remains the same.

It's worth noting that binning can introduce errors in the spectral data,
especially if the signal-to-noise ratio is low or if the binning is too
coarse. In general, astronomers choose a binning size that balances the need
for a high signal-to-noise ratio with the desire to maintain the spectral
resolution and avoid introducing significant errors in the data.

In summary, the principle of flux density conservation is important to
consider when binning astronomical spectra, and astronomers need to scale the
intensity of each binned pixel to ensure that the total energy emitted by the
object is conserved.

Here, the method used is with the Cumulative function to produce a new
flux-conserved, some options can be chosen for the interpolation:

<center>
<table>

<tr><td>Integer Number</td><td>Option: Interpolation Schemes</td><td>Brief Description</td></tr>

<tr><td>0)</td><td>LINinterpol</td><td>BBBBBBBBB</td></tr>
<tr><td>1)</td><td>SPLINE3DFor</td><td>BBBBBBBBB</td></tr>
<tr><td>2)</td><td>SPLINE1DArr</td><td>BBBBBBBBB</td></tr>
<tr><td>3)</td><td>AkimaSpline</td><td>BBBBBBBBB</td></tr>
<tr><td>4)</td><td>Interpolado</td><td>BBBBBBBBB</td></tr>
<tr><td>5)</td><td>LINdexerpol</td><td>BBBBBBBBB</td></tr>

</table>
</center>

You can easily install <a
href=https://pypi.org/project/pyfluxconserving/>pyfluxconserving</a> by using pip -
<a href='https://pypi.org/'>PyPI - The Python Package Index</a>: <pre> <code>
pip install pyfluxconserving </code> </pre> or by using a generated conda
repository <a
href='https://anaconda.org/neutrinomuon/pyfluxconserving'>https://anaconda.org/neutrinomuon/pyfluxconserving</a>:

[![badgetanaconda](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/version.svg)](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/version.svg)
[![badgetreleasedate](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/latest_release_date.svg)](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/latest_release_date.svg)
[![badgetplatforms](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/platforms.svg
)](https://anaconda.org/neutrinomuon/pyfluxconserving/badges/platforms.svg)
<pre>
<code>
conda install -c neutrinomuon pyfluxconserving
</code>
</pre>
OBS.: Linux, OS-X and Windows pre-compilations available in conda.

You can also clone the repository and install by yourself in your machine:
<pre>
<code>
git clone https://github.com/neutrinomuon/FluxConserving
python setup.py install
</code>
</pre>

The methods are given by Int_Type and may be summarized bellow:

<table>
<tr><td>Int_Type</td><td>Type</td><td>Description</td></tr>
<tr><td>0<td>R</td><td>Right rectangle Integral  </td></tr>
<tr><td>1<td>L</td><td>Left rectangle Integral   </td></tr>
<tr><td>2<td>T</td><td>Trapezoidal rule          </td></tr>
<tr><td>3<td>S</td><td>Simple Integral           </td></tr>
<tr><td>4<td>M</td><td>Median rectangle Integral </td></tr>
<tr><td>5<td>I</td><td>Simpsonregel's rule       </td></tr>
<tr><td>6<td>G</td><td>Gauss-Legendre Quadrature </td></tr>
</table>

The main structure of the directories and files are:

<pre>
<code>
IntegralALL
├── dist
│   └── pyintegralall-0.0.1.tar.gz
├── README.md
├── pyintegralall.egg-info
│   ├── PKG-INFO
│   ├── dependency_links.txt
│   ├── SOURCES.txt
│   ├── top_level.txt
│   └── requires.txt
├── LICENSE.txt
├── setup.py
├── build.bat
├── tutorials
│   ├── Example1 - IntegralALL.ipynb
│   └── .ipynb_checkpoints
│       └── Example1 - IntegralALL-checkpoint.ipynb
├── build.sh
├── src
│   ├── python
│   │   ├── __init__.py
│   │   └── PyIntegralALL.py
│   └── fortran
│       ├── IntegralALL.compile
│       ├── IntegralALL.f90
│       ├── IntegralALL.cpython-39-x86_64-linux-gnu.so
│       ├── DataTypes.f90
│       ├── GaussLegendreQuadrature.cpython-39-x86_64-linux-gnu.so
│       ├── LINinterpol.compile
│       ├── GaussLegendreQuadrature.f90
│       ├── LINinterpol.cpython-39-x86_64-linux-gnu.so
│       ├── LINinterpol.f90
│       └── GaussLegendreQuadrature.compile
├── version.txt
├── meta.yaml
└── build
    ├── lib.linux-x86_64-3.9
    │   └── pyintegralall
    ├── src.linux-x86_64-3.9
    │   ├── pyintegralall
    │   ├── build
    │   └── numpy
    └── temp.linux-x86_64-3.9
        ├── pyintegralall
        ├── __pycache__
        ├── ccompiler_opt_cache_ext.py
        ├── src
        └── build

19 directories, 28 files
</code>
</pre>

PyFluxConSPec.py is a python wrapper to the library in fortran called
pyfluxconserving.flib. The fortran directory can be compiled separately for
each individual subroutine.
