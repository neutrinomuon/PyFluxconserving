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

<tr><td>0)</td><td>AkimaSpline</td><td>Akima Spline interpolation. The Akima
spline is a C1 differentiable function (that is, has a continuous first
derivative) but, in general, will have a discontinuous second derivative at
the knot points.<p>
                                                                           
<ol>
    <li>Akima, H. (1970). A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures. <i>Journal of the ACM</i>, 17(4), 589-602. Association for Computing Machinery. DOI: 10.1145/321607.321609. Link: <a href="https://dl.acm.org/doi/10.1145/321607.321609">https://dl.acm.org/doi/10.1145/321607.321609</a></li>
    <li>De Boor, C. (1978). A Practical Guide to Splines. Springer-Verlag. DOI: 10.1007/978-1-4757-0457-4. Link: <a href="https://www.springer.com/gp/book/9780387953663">https://www.springer.com/gp/book/9780387953663</a></li>
    <li>Forsythe, G.E. (1979) <i>Computer Methods for Mathematical Computations</i>. Prentice-Hall, Inc. Link: <a href="https://onlinelibrary.wiley.com/doi/10.1002/zamm.19790590235">https://onlinelibrary.wiley.com/doi/10.1002/zamm.19790590235</a></li>
    <li>Bartels, R.H., Beatty, J.C., and Barsky, B.A. (1998). An Introduction to Splines for Use in Computer Graphics and Geometric Modeling. Morgan Kaufmann Publishers. ISBN: 978-1558604000. Link: <a href="https://www.sciencedirect.com/science/article/abs/pii/0378475494001228">https://www.sciencedirect.com/science/article/abs/pii/0378475494001228</a></li>
    <li>Cheney, W. and Kincaid, D. (2008). Numerical Mathematics and Computing (6th ed.). Brooks/Cole. ISBN: 978-0495114758. Link: <a href="https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758">https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758</a></li>
    <li>Burden, R. L., &amp; Faires, J. D. (2010). Numerical Analysis (9th ed.). Brooks/Cole. ISBN: 978-0538733519. Link: <a href="https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/0538733519">https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/0538733519</a></li>
    <li>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery, B. P. (2007). Numerical Recipes: The Art of Scientific Computing (3rd ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a href="http://numerical.recipes/">http://numerical.recipes/</a></li>
</ol>

</td></tr>

<tr><td>1)</td><td>Interpolado</td><td>Based on a linear interpolation within
a table of pair values.<p>

<ol>

<li>Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery,
B. P. (2007). <i>Numerical recipes: The art of scientific computing</i> (3rd
ed.). Cambridge University Press. DOI: <a
href="https://doi.org/10.1017/CBO9780511807447">10.1017/CBO9780511807447</a></li>

<li>Burden, R. L., & Faires, J. D. (2010). <i>Numerical analysis</i> (9th
ed.). Brooks/Cole. DOI: <a
href="https://doi.org/10.1007/978-0-387-77854-3">10.1007/978-0-387-77854-3</a></li>

<li>Chapra, S. C., & Canale, R. P. (2010). <i>Numerical methods for
engineers</i> (6th ed.). McGraw-Hill. DOI: <a
href="https://doi.org/10.1080/03081079808818664">10.1080/03081079808818664</a></li>

<li>Atkinson, K. E. (1989). <i>An introduction to numerical analysis</i> (2nd
ed.). John Wiley & Sons. DOI: <a
href="https://doi.org/10.1002/9780470753763">10.1002/9780470753763</a></li>

<li>Burden, R. L., & Faires, J. D. (1985). <i>Numerical analysis</i> (3rd
ed.). PWS-Kent. DOI: <a
href="https://doi.org/10.2307/2321069">10.2307/2321069</a></li>
</ol>

</td></tr>

<tr><td>2)</td><td>LINdexerpol</td><td>Based on a linear interpolation within
a table of pair values using indexing. The references are the same as in
1).</td></tr>

<tr><td>3)</td><td>LINinterpol</td><td>Based on a linear interpolation within
a table of pair values. The references are the same as in 1).</td></tr>

<tr><td>3)</td><td>LINinterpol</td><td>Based on a linear interpolation within
a table of pair values. The references are the same as in 1).</td></tr>

<tr><td>4)</td><td>SPLINE1DArr</td><td>This is a Fortran 2003 subroutine
called SPLINE1DArr that takes an array of values x to interpolate from the
arrays t and y. It has ten input arguments, six output arguments, and two
optional arguments. The interpolation is linear.<p>

<ol>

<li>Forsythe, G.E. (1977) <em>Computer Methods For Mathematical
Computations</em>. Ed. Prentice-Hall, Inc. <a
href="https://www.abebooks.com/9780131653323/Computer-Methods-Mathematical-Computations-Prentice-Hall-0131653328/plp">Link</a></li>

<li>De Boor, C. (1978). <em>A Practical Guide to
Splines</em>. Springer-Verlag. DOI: <a
href="https://doi.org/10.1007/978-1-4757-0457-4">10.1007/978-1-4757-0457-4</a>. Link:
<a href="https://www.springer.com/gp/book/9780387953663">Link</a></li>

<li>Smith, L. (1995). A cubic spline algorithm based on a modified divided
difference. <em>Mathematics and Computers in Simulation</em>, 39(3),
347-350. DOI: <a
href="https://doi.org/10.1016/0378-4754(94)00122-8">10.1016/0378-4754(94)00122-8</a>. Link:
<a
href="https://www.sciencedirect.com/science/article/abs/pii/0378475494001228">Link</a></li>

<li>Bartels, R.H., Beatty, J.C., and Barsky, B.A. (1998). <em>An Introduction
to Splines for Use in Computer Graphics and Geometric Modeling</em>. Morgan
Kaufmann Publishers. ISBN: 978-1558604000. Link: <a
href="https://www.elsevier.com/books/an-introduction-to-splines-for-use-in-computer-graphics-and-geometric-modeling/bartels/978-1-55860-400-0">Link</a></li>

<li>Burden, R. L., &amp; Faires, J. D. (2010). <em>Numerical Analysis</em>
(9th ed.). Brooks/Cole. ISBN: 978-0538733519. Link: <a
href="https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/0538733519">Link</a></li>

<li>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery,
B. P. (2007). <em>Numerical Recipes: The Art of Scientific Computing</em> (3rd
ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a
href="http://numerical.recipes/">Link</a></li>

<li>Cheney, W. and Kincaid, D. (2008). <em>Numerical Mathematics and
Computing</em> (6th ed.). Brooks/Cole. ISBN: 978-0495114758. Link: <a
href="https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758">Link</a></li>

</ol>
</td></tr>

<tr><td>5)</td><td>SPLINE3DFor</td><td>This function evaluates the cubic
spline interpolation. The same references as in 4)</td></tr>

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

The main structure of the directories and files are:

<pre>
<code>
FluxConserving
├── dist
│   └── pyfluxconserving-0.0.1.tar.gz
├── README.md
├── LICENSE.txt
├── setup.py
├── tutorials
│   ├── fluxconserving.png
│   └── Flux-Conserving Example.ipynb
├── pyfluxconserving.egg-info
│   ├── PKG-INFO
│   ├── dependency_links.txt
│   ├── SOURCES.txt
│   ├── top_level.txt
│   └── requires.txt
├── src
│   ├── python
│   │   ├── PyFluxConSpec.py
│   │   ├── califa_cmap_alternative.py
│   │   ├── PyLinear__int.py
│   │   ├── PyLINinterpol.py
│   │   ├── PySPLINECubic.py
│   │   ├── fluxconserve.py
│   │   ├── __init__.py
│   │   ├── PySPLINE1DArr.py
│   │   ├── PySPLINE3DFor.py
│   │   ├── fluxconserving.png
│   │   ├── PyAkimaSpline.py
│   │   ├── PySPLINE3DAlt.py
│   │   ├── PyInterpolado.py
│   │   └── PyLINdexerpol.py
│   └── fortran
│       ├── LINdexerpol.f90
│       ├── SPLINE3DFor.cpython-39-x86_64-linux-gnu.so
│       ├── SPLINE1DFlt.cpython-39-x86_64-linux-gnu.so
│       ├── SPLINE3DFor.compile
│       ├── FluxConSpec.compile
│       ├── SPLINE1DArr.compile
│       ├── Interpolado.cpython-39-x86_64-linux-gnu.so
│       ├── Interpolado.cpython-38-x86_64-linux-gnu.so
│       ├── SPLINE1DFlt.f90
│       ├── AkimaSpline.f90
│       ├── SPLINE1DArr.f90
│       ├── SPLINE1DFlt.cpython-38-x86_64-linux-gnu.so
│       ├── SPLINE3DFor.cpython-38-x86_64-linux-gnu.so
│       ├── AkimaSpline.compile
│       ├── DataTypes.f90
│       ├── LINdexerpol.compile
│       ├── SPLINE3DFor.f90
│       ├── SPLINE1DFlt.compile
│       ├── Interpolado.f90
│       ├── SPLINE1DFlt.cpython-310-x86_64-linux-gnu.so
│       ├── FluxConSpec.f90
│       ├── Interpolado.cpython-310-x86_64-linux-gnu.so
│       ├── LINdexerpol.cpython-310-x86_64-linux-gnu.so
│       ├── README.txt
│       ├── Interpolado.compile
│       ├── SPLINE3DFor.cpython-310-x86_64-linux-gnu.so
│       ├── LINinterpol.f90
│       ├── LINdexerpol.cpython-39-x86_64-linux-gnu.so
│       └── LINdexerpol.cpython-38-x86_64-linux-gnu.so
├── version.txt
└── build
    ├── lib.linux-x86_64-3.9
    │   └── pyfluxconserving
    ├── src.linux-x86_64-3.9
    │   ├── pyfluxconserving
    │   ├── build
    │   └── numpy
    └── temp.linux-x86_64-3.9
        ├── pyfluxconserving
        ├── ccompiler_opt_cache_ext.py
        ├── src
        ├── .libs
        └── build

18 directories, 56 files
</code>
</pre>

PyFluxConSPec.py is a python wrapper to the library in fortran called
pyfluxconserving.flib. The fortran directory can be compiled separately for
each individual subroutine.
