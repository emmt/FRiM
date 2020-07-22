# FRactal Iterative Method

FRiM (FRactal Iterative Method) is an algorithm to generate/recover random
fields with stationnary covariance like turbulent wavefronts.

This repository provides the code corresponding to the method described in
Thiébaut & Tallon, *"Fast minimum variance wavefront reconstruction for
extremely large telescopes"*, J. Opt. Soc. Am. A, vol. **27**, pp. 1046-1059
(2010).

This repository is organized as follows:
* [`c`](./c) contains the C sources;
* [`matlab-octave`](./matlab-octave) contains code to generate turbulent
  wavefronts with [MATLAB](https://en.wikipedia.org/wiki/MATLAB) or
  [Octave](https://www.gnu.org/software/octave/);
* [`yorick`](./yorick) contains the source of the FRiM plugin for
  [Yorick](http://github.com/dhmunro/yorick).


## Yorick plug-in

To build and install the [Yorick](http://github.com/dhmunro/yorick) plug-in:

```.sh
cd $BUILD
$SRCDIR/configure
make clean
make
make install
```

where `$BUILD` is the directory where to build the plug-in and
`$SRCDIR/configure` is the path to the configuration script in the
[`yorick`](./yorick) directory.  The directories `$BUILD` and `$SRCDIR` can be
the same.
