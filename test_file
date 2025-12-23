### PyFluxconserving
####  Fully based on a Fortran legacy package to easily compute the flux-density conservation
Obs.: A Fortran legacy Interpolation routines also furnished in this package <br>
email: [antineutrinomuon@gmail.com](mailto:antineutrinomuon@gmail.com)

linkedin: <a href="https://www.linkedin.com/in/jean-michel-gomes/">https://www.linkedin.com/in/jean-michel-gomes/</a>

github repository: <a href="https://github.com/neutrinomuon/PyFluxconserving">PyFluxconserving</a>

last stable version: 0.0.16

© Copyright ®

J.G. - Jean Gomes @ 2025

<hr>

<img src="https://skillicons.dev/icons?i=python,fortran,c,numpy&theme=light"><br>
<img src="https://img.shields.io/pypi/pyversions/PyFluxconserving"><img src="https://anaconda.org/neutrinomuon/PyFluxconserving/badges/license.svg">

<hr>

<div align="center">
<img src='https://raw.githubusercontent.com/neutrinomuon/PyFluxconserving/main/figures/PyFluxconserving_.png' width='60%'>
</div>

<hr>

#### <b>RESUME</b>

<img
src="https://raw.githubusercontent.com/neutrinomuon/PyFluxconserving/main/figures/PyFluxconserving.png"
width=120>

Original Fortran 2003+ legacy routines date back to 2003–2004. Please
refer to the <a
href='https://github.com/neutrinomuon/PyFluxconserving/blob/main/LICENSE.txt'>LICENSE.txt</a>
file. These routines have been converted into Python libraries.

When analysing astronomical spectra, astronomers often bin the data to
increase the signal-to-noise ratio and mitigate the effects of
noise. Binning refers to the process of averaging the intensity of
adjacent spectral channels, or pixels, to produce a new, coarser
dataset.

During binning, it is crucial to maintain the principle of flux
density conservation. This means that the total energy emitted by the
object, and consequently its flux density, must remain constant after
binning.

To conserve flux density, the intensity of each binned pixel should be
scaled according to the number of pixels it represents. For instance,
if two adjacent pixels are binned together, the intensity of the
resulting bin should equal the sum of the intensities of the two
original pixels, divided by two. This ensures that the total energy in
the bin is conserved and that the object's flux density remains
unchanged.

It is important to note that binning can introduce errors in spectral
data, especially if the signal-to-noise ratio is low or if the binning
is too coarse. Nevertheless, this procedure is often necessary, as the
observed spectrum of an astronomical source is generally
non-linear. Astronomers usually choose a binning size that balances
the need for a high signal-to-noise ratio with the requirement to
maintain spectral resolution and minimise significant errors in the
data.

In summary, flux density conservation is a key consideration when
binning astronomical spectra. Astronomers must scale the intensity of
each binned pixel to ensure that the total energy emitted by the
object is preserved.

Note: SpectRes by A. Carnall is not included in this distribution but
is provided for comparison: <a
href='https://github.com/ACCarnall/SpectRes'>https://github.com/ACCarnall/SpectRes
</a>. To install it for comparison, please refer to the repository.
However, it is not necessary for the usage of this
package. Accompanying there are several routines for interpolations.

<hr>

#### <b>INSTALLATION</b>

You can easily install <a
href=https://pypi.org/project/PyFluxconserving/>PyFluxconserving</a> by using pip -
<a href='https://pypi.org/'>PyPI - The Python Package Index</a>:

<pre>
pip install PyFluxconserving
</pre>

<br>or by using a generated conda repository <a
href='https://anaconda.org/neutrinomuon/PyFluxconserving'>https://anaconda.org/neutrinomuon/PyFluxconserving</a>:

<img src="https://anaconda.org/neutrinomuon/PyFluxconserving/badges/version.svg"><img src="https://anaconda.org/neutrinomuon/PyFluxconserving/badges/latest_release_date.svg"><img src="https://anaconda.org/neutrinomuon/PyFluxconserving/badges/platforms.svg">

<pre>
conda install -c neutrinomuon PyFluxconserving
</pre>

<br>OBS.: Linux, OS-X and Windows pre-compilations available in conda.

You can also clone the repository and install by yourself in your machine:

<pre>
git clone https://github.com/neutrinomuon/PyFluxconserving
pip install --user ninja (optional)
pip install --user meson (optional)
rm -rf build *.egg-info PyFluxconserving
python setup.py build_ext --inplace
pip install -e .
</pre>

<hr>

#### <b>METHOD & REFERENCES</b>

Here, the method used is with the Cumulative function to produce a new
flux-conserved, some options can be chosen for the interpolation:

<table border="1">

<tr><td width='20%'>Integer Number</td><td width='33%'>Option: Interpolation Schemes</td><td>Brief Description</td></tr>

<tr><td>0)</td><td>AkimaSpline</td><td>Akima Spline interpolation. The Akima spline is a C1 differentiable function (that is, has a continuous first
derivative) but, in general, will have a discontinuous second derivative at the knot points.<p>

<ul>
<il>Akima, H. (1970). A New Method of Interpolation and Smooth Curve Fitting Based on Local Procedures. <i>Journal of the ACM</i>, 17(4),
589-602. Association for Computing Machinery. DOI: 10.1145/321607.321609. Link: <a href="https://dl.acm.org/doi/10.1145/321607.321609">https://dl.acm.org/doi/10.1145/321607.321609</a></il><p>

<il>De Boor, C. (1978). A Practical Guide to Splines. Springer-Verlag. ISBN: 978-0-387-95366-3. Link: <a href="https://www.springer.com/gp/book/9780387953663">https://www.springer.com/gp/book/9780387953663</a></il><p>

<il>Forsythe, G.E. (1979) <i>Computer Methods for Mathematical Computations</i>. Prentice-Hall, Inc. DOI: 10.1002/zamm.19790590235. Link: <a href="https://onlinelibrary.wiley.com/doi/10.1002/zamm.19790590235">https://onlinelibrary.wiley.com/doi/10.1002/zamm.19790590235</a></il><p>

<il>Bartels, R.H., Beatty, J.C., and Barsky, B.A. (1987). An Introduction to Splines for Use in Computer Graphics and Geometric Modeling. Morgan Kaufmann
Publishers. Link: <a href="https://www.osti.gov/biblio/5545263">https://www.osti.gov/biblio/5545263</a></il><p>

<il>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery, B. P. (2007). Numerical Recipes: The Art of Scientific Computing (3rd
ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a href="http://numerical.recipes/">http://numerical.recipes/</a></il><p>

<il>Cheney, W. and Kincaid, D. (2008). Numerical Mathematics and Computing (6th ed.). Brooks/Cole. ISBN: 978-0495114758. Link: <a href="https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758">https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758</a></il><p>

<il>Burden, R. L., &amp; Faires, J. D. (2010). Numerical Analysis (9th ed.). Brooks/Cole. ISBN: 978-0538733519. Link: <a href="https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/0538733519">https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/0538733519</a></il><p>

</ul>

</td></tr>

<tr><td>1)</td><td>Interpolado</td><td>Based on a linear interpolation within
a table of pair values.<p>

<ul>

<il>Atkinson, K. E. (1991). <i>An introduction to numerical analysis</i> (2nd ed.). John Wiley & Sons. ISBN: 978-0-471-62489-9. Link: <a href="https://www.wiley.com/en-us/An+Introduction+to+Numerical+Analysis%2C+2nd+Edition-p-9780471624899">https://www.wiley.com/en-us/An+Introduction+to+Numerical+Analysis%2C+2nd+Edition-p-9780471624899</a></il><p>

<il>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery, B. P. (2007). <em>Numerical Recipes: The Art of Scientific Computing</em> (3rd
ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a href="http://numerical.recipes/">http://numerical.recipes/</a></il><p>

<il>Chapra, S. C., & Canale, R. P. (2010). <i>Numerical methods for engineers</i> (6th ed.). McGraw-Hill. ISBN: 978-0073401065. Link: <a href="https://www.amazon.com/Numerical-Methods-Engineers-Steven-Chapra/dp/0073401064">https://www.amazon.com/Numerical-Methods-Engineers-Steven-Chapra/dp/0073401064</a></il><p>

<il>Burden, R. L., Faires, J. D. & A.M. Burden (2015). <i>Numerical analysis (10th ed.)</i> (9th ed.). Cengage Learning. ISBN: 978-1305253667. Link: <a href="https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/1305253663">https://www.amazon.com/Numerical-Analysis-Richard-L-Burden/dp/1305253663</a></il><p>

</ul>

</td></tr>

<tr><td>2)</td><td>LINdexerpol</td><td>Based on a linear interpolation within
a table of pair values using indexing. The references are the same as in
1).</td></tr>

<tr><td>3)</td><td>LINinterpol</td><td>Based on a linear interpolation within
a table of pair values. The references are the same as in 1).</td></tr>

<tr><td>4)</td><td>SPLINE1DArr</td><td>This is a Fortran 2003 subroutine
called SPLINE1DArr that takes an array of values x to interpolate from the
arrays t and y. It has ten input arguments, six output arguments, and two
optional arguments. The interpolation is linear.<p>

<ul>
<il>Forsythe, G.E. (1977) <em>Computer Methods For Mathematical Computations</em>. Ed. Prentice-Hall, Inc. <a href="https://onlinelibrary.wiley.com/doi/abs/10.1002/zamm.19790590235">https://onlinelibrary.wiley.com/doi/abs/10.1002/zamm.19790590235</a></il><p>

<il>De Boor, C. (1978). <em>A Practical Guide to Splines</em>. Springer-Verlag. ISBN: 978-0-387-95366-3. Link: <a href="https://www.springer.com/gp/book/9780387953663">https://www.springer.com/gp/book/9780387953663</a></il><p>

<il>Bartels, R.H., Beatty, J.C., and Barsky, B.A. (1998). <em>An Introductionto Splines for Use in Computer Graphics and Geometric Modeling</em>. Morgan
Kaufmann Publishers. ISBN: 978-1558604000. Link: <a href="https://www.amazon.com/Introduction-Computer-Graphics-Geometric-Modeling/dp/1558604006">https://www.amazon.com/Introduction-Computer-Graphics-Geometric-Modeling/dp/1558604006</a></il><p>

<il>Press, W. H., Teukolsky, S. A., Vetterling, W. T., &amp; Flannery, B. P. (2007). <em>Numerical Recipes: The Art of Scientific Computing</em> (3rd
ed.). Cambridge University Press. ISBN: 978-0521880688. Link: <a href="http://numerical.recipes/">http://numerical.recipes/</a></il><p>

<il>Cheney, W. and Kincaid, D. (2007). <em>Numerical Mathematics and Computing</em> (6th ed.). Brooks/Cole. ISBN: 978-0495114758. Link: <a href="https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758">https://www.amazon.com/Numerical-Mathematics-Computing-Ward-Cheney/dp/0495114758</a></il><p>

</td></tr>

<tr><td>5)</td><td>SPLINE3DFor</td><td>This function evaluates the cubic
spline interpolation. The same references as in 4).</td></tr>

</table>

<hr>

#### <b>STRUCTURE</b>
<pre>
#################################################
workspace
├── version.txt
├── .github
│   └── workflows
│       ├── main.yml
│       └── pylint.yml
├── src
│   ├── fortran
│   │   ├── Interpolado.compile
│   │   ├── SPLINE3DFor.cpython-310-x86_64-linux-gnu.so
│   │   ├── LINdexerpol.cpython-311-darwin.so
│   │   ├── SPLINE3DFor.cpython-310-darwin.so
│   │   ├── LINinterpol.compile
│   │   ├── AkimaSpline.cpython-311-darwin.so
│   │   ├── Interpolado.f90
│   │   ├── FluxConSpec.cpython-311-darwin.so
│   │   ├── Interpolado.cpython-310-x86_64-linux-gnu.so
│   │   ├── SPLINE1DArr.f90
│   │   ├── SPLINE1DFlt.f90
│   │   ├── Interpolado.cpython-39-x86_64-linux-gnu.so
│   │   ├── FluxConSpec.cpython-310-x86_64-linux-gnu.so
│   │   ├── SPLINE1DArr.cpython-39-darwin.so
│   │   ├── SPLINE1DFlt.cpython-310-x86_64-linux-gnu.so
│   │   ├── AkimaSpline.cpython-312-x86_64-linux-gnu.so
│   │   ├── SPLINE1DArr.cpython-310-darwin.so
│   │   ├── LINdexerpol.cpython-310-x86_64-linux-gnu.so
│   │   ├── AkimaSpline.cpython-310-darwin.so
│   │   ├── Interpolado.cpython-38-x86_64-linux-gnu.so
│   │   ├── SPLINE3DFor.f90
│   │   ├── LINinterpol.cpython-311-x86_64-linux-gnu.so
│   │   ├── SPLINE3DFor.cpython-39-darwin.so
│   │   ├── SPLINE1DFlt.compile
│   │   ├── FluxConSpec.f90
│   │   ├── FluxConSpec.cpython-39-x86_64-linux-gnu.so
│   │   ├── FluxConSpec.cpython-310-darwin.so
│   │   ├── SPLINE3DFor.cpython-311-darwin.so
│   │   ├── LINinterpol.f90
│   │   ├── SPLINE1DFlt.cpython-39-x86_64-linux-gnu.so
│   │   ├── AkimaSpline.compile
│   │   ├── SPLINE3DFor.compile
│   │   ├── AkimaSpline.cpython-39-darwin.so
│   │   ├── LINinterpol.cpython-310-darwin.so
│   │   ├── FluxConSpec.cpython-39-darwin.so
│   │   ├── LINinterpol.cpython-310-x86_64-linux-gnu.so
│   │   ├── LINdexerpol.f90
│   │   ├── SPLINE1DFlt.cpython-311-darwin.so
│   │   ├── SPLINE1DFlt.cpython-39-darwin.so
│   │   ├── README.txt
│   │   ├── SPLINE3DFor.cpython-38-x86_64-linux-gnu.so
│   │   ├── FluxConSpec.compile
│   │   ├── DataTypes.f90
│   │   ├── LINdexerpol.compile
│   │   ├── Interpolado.cpython-312-x86_64-linux-gnu.so
│   │   ├── LINdexerpol.cpython-38-x86_64-linux-gnu.so
│   │   ├── SPLINE1DArr.compile
│   │   ├── LINdexerpol.cpython-310-darwin.so
│   │   ├── SPLINE1DFlt.cpython-38-x86_64-linux-gnu.so
│   │   ├── LINinterpol.cpython-311-darwin.so
│   │   ├── AkimaSpline.cpython-38-x86_64-linux-gnu.so
│   │   ├── FluxConSpec.cpython-38-x86_64-linux-gnu.so
│   │   ├── SPLINE3DFor.cpython-39-x86_64-linux-gnu.so
│   │   ├── SPLINE1DArr.cpython-311-darwin.so
│   │   ├── AkimaSpline.f90
│   │   ├── SPLINE1DFlt.cpython-310-darwin.so
│   │   ├── LINdexerpol.cpython-39-darwin.so
│   │   ├── LINinterpol.cpython-39-darwin.so
│   │   └── LINdexerpol.cpython-39-x86_64-linux-gnu.so
│   └── python
│       ├── __init__.py
│       ├── PyAkimaSpline.py
│       ├── fluxconserve.py
│       ├── PyInterpolado.py
│       ├── PySPLINE1DArr.py
│       ├── PyLinear__int.py
│       ├── PyLINdexerpol.py
│       ├── califa_cmap_alternative.py
│       ├── PySPLINE3DFor.py
│       ├── fluxconserving.png
│       ├── PyFluxConSpec.py
│       └── PyLINinterpol.py
├── scripts
│   └── update_readme.py
├── build
│   ├── lib.macosx-11.0-arm64-cpython-39
│   │   └── PyFluxconserving
│   │       ├── __init__.py
│   │       ├── PyAkimaSpline.py
│   │       ├── fluxconserve.py
│   │       ├── PyInterpolado.py
│   │       ├── PySPLINE1DArr.py
│   │       ├── PyLinear__int.py
│   │       ├── PyLINdexerpol.py
│   │       ├── flib.cpython-39-darwin.so
│   │       ├── califa_cmap_alternative.py
│   │       ├── PySPLINE3DFor.py
│   │       ├── PyFluxConSpec.py
│   │       └── PyLINinterpol.py
│   ├── src.macosx-11.0-arm64-3.9
│   │   ├── build
│   │   │   └── src.macosx-11.0-arm64-3.9
│   │   │       └── PyFluxconserving
│   │   │           ├── fortranobject.c
│   │   │           └── fortranobject.h
│   │   ├── numpy
│   │   │   └── distutils
│   │   │       └── include
│   │   │           └── npy_cpu_dispatch_config.h
│   │   └── PyFluxconserving
│   │       ├── flib-f2pywrappers.f
│   │       ├── flibmodule.c
│   │       └── flib-f2pywrappers2.f90
│   └── temp.macosx-11.0-arm64-cpython-39
│       ├── src
│       │   └── fortran
│       │       ├── LINdexerpol.o
│       │       ├── DataTypes.o
│       │       ├── LINinterpol.o
│       │       ├── FluxConSpec.o
│       │       ├── AkimaSpline.o
│       │       ├── SPLINE1DArr.o
│       │       ├── Interpolado.o
│       │       └── SPLINE3DFor.o
│       ├── build
│       │   └── src.macosx-11.0-arm64-3.9
│       │       ├── build
│       │       │   └── src.macosx-11.0-arm64-3.9
│       │       │       └── PyFluxconserving
│       │       │           ├── fortranobject.o
│       │       │           └── fortranobject.o.d
│       │       └── PyFluxconserving
│       │           ├── flib-f2pywrappers2.o
│       │           ├── flibmodule.o
│       │           ├── flibmodule.o.d
│       │           └── flib-f2pywrappers.o
│       ├── __pycache__
│       │   └── ccompiler_opt_cache_ext.cpython-39.pyc
│       ├── ccompiler_opt_cache_ext.py
│       └── PyFluxconserving
│           └── moddatatype.mod
├── dist
│   ├── pyfluxconserving-0.1.16-py3-none-any.whl
│   └── pyfluxconserving-0.1.16.tar.gz
├── requirements.txt
├── requirements.py
├── README.md
├── setup.py
├── tutorials
│   ├── PyFluxconserving Example.ipynb
│   ├── .virtual_documents
│   │   └── PyFluxconserving Example.ipynb
│   └── .ipynb_checkpoints
│       ├── pyfluxconserving Example-checkpoint.ipynb
│       └── PyFluxconserving Example-checkpoint.ipynb
├── figures
│   ├── PyFluxconserving_.png
│   └── PyFluxconserving.png
├── showdown.min.js
├── index.html
├── LICENSE.txt
├── PyFluxconserving
│   ├── __init__.py
│   ├── PyAkimaSpline.py
│   ├── fluxconserve.py
│   ├── PyInterpolado.py
│   ├── PySPLINE1DArr.py
│   ├── PyLinear__int.py
│   ├── __pycache__
│   │   ├── PyFluxConSpec.cpython-312.pyc
│   │   ├── califa_cmap_alternative.cpython-312.pyc
│   │   ├── __init__.cpython-312.pyc
│   │   └── fluxconserve.cpython-312.pyc
│   ├── flib.pyf
│   ├── flib.cpython-312-x86_64-linux-gnu.so
│   ├── PyLINdexerpol.py
│   ├── califa_cmap_alternative.py
│   ├── PySPLINE3DFor.py
│   ├── PyFluxConSpec.py
│   └── PyLINinterpol.py
├── README_setup.txt
├── .git
│   ├── logs
│   │   ├── HEAD
│   │   └── refs
│   │       ├── remotes
│   │       │   └── origin
│   │       │       └── main
│   │       └── heads
│   │           └── main
│   ├── config
│   ├── HEAD
│   ├── refs
│   │   ├── remotes
│   │   │   └── origin
│   │   │       └── main
│   │   ├── heads
│   │   │   └── main
│   │   └── tags
│   ├── config.worktree
│   ├── FETCH_HEAD
│   ├── description
│   ├── hooks
│   │   ├── sendemail-validate.sample
│   │   ├── pre-receive.sample
│   │   ├── fsmonitor-watchman.sample
│   │   ├── prepare-commit-msg.sample
│   │   ├── post-update.sample
│   │   ├── pre-rebase.sample
│   │   ├── push-to-checkout.sample
│   │   ├── pre-merge-commit.sample
│   │   ├── pre-applypatch.sample
│   │   ├── pre-commit.sample
│   │   ├── commit-msg.sample
│   │   ├── applypatch-msg.sample
│   │   ├── pre-push.sample
│   │   └── update.sample
│   ├── index
│   ├── objects
│   │   ├── pack
│   │   │   ├── pack-46f3154c766289426b0cda7512267f9284038e5e.rev
│   │   │   ├── pack-46f3154c766289426b0cda7512267f9284038e5e.idx
│   │   │   └── pack-46f3154c766289426b0cda7512267f9284038e5e.pack
│   │   └── info
│   ├── shallow
│   └── info
│       └── exclude
├── pyfluxconserving
│   └── meta.yaml
├── tree.out
└── PyFluxconserving.egg-info
    ├── PKG-INFO
    ├── dependency_links.txt
    ├── top_level.txt
    └── SOURCES.txt

53 directories (0 symlink), 179 files (0 symlink)
#################################################
Generated: treehue_colored @2025 - © Jean Gomes -
#################################################
</pre>

<br>PyFluxConSPec.py is a python wrapper to the library in fortran called
PyFluxconserving.flib. The fortran directory can be compiled separately for
each individual subroutine.

<hr>

#### ISSUES AND CONTRIBUTIONS

If you encounter any issues with this project, please feel free to submit an issue on the GitHub repository. We appreciate your feedback and are committed to improving the quality of our codebase.

If you'd like to contribute to this project, we welcome pull requests from the community. Before submitting a pull request, please make sure to fork the repository and create a new branch for your changes. Once your changes are complete, submit a pull request and we'll review your code as soon as possible.

For any questions or concerns about contributing, please contact the project maintainer at antineutrinomuon@gmail.com. Thank you for your interest in contributing to our project!

<hr>

#### <b>LICENSE</b>

This software is provided "AS IS" (see DISCLAIMER below). Permission to
use, for non-commercial purposes is granted. Permission to modify for personal
or internal use is granted, provided this copyright and disclaimer are
included in ALL copies of the software. All other rights are reserved. In
particular, redistribution of the code is not allowed without explicit
permission by the author.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
