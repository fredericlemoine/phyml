
[![Tests](https://github.com/fredericlemoine/phyml/actions/workflows/tests.yaml/badge.svg?branch=dev_fred_struct)](https://github.com/fredericlemoine/phyml/actions/workflows/tests.yaml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/phyml/README.html) 


### Overview

PhyML is a software package that uses modern statistical approaches to analyse alignments of nucleotide or amino acid sequences in a phylogenetic
framework. The main tool in this package builds phylogenies under the maximum likelihood criterion. It implements a large number of substitution
models coupled to efficient options to search the space of phylogenetic tree topologies. PhyTime is another tool in the PhyML package that focuses
on divergence date estimation in a Bayesian setting. The main strengths of PhyTime lies in its ability to accommodate for uncertrainty in the placement of fossil calibration
and the use of realistic models of rate variation along the tree. Finally, PhyREX fits the spatial-Lambda-Fleming-Viot
model to geo-referenced genetic data. This model is similar to the structured coalescent, but assumes that individuals are distributed along a spatial continuum rather
than discrete demes. PhyREX can be used to estimate population densities and rates of dispersal. Its output can be processed by treeannotator (from the BEAST package) as well as
SPREAD.

### Citations

- New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0 S Guindon, JF Dufayard, V Lefort, M Anisimova, W Hordijk, O Gascuel *Systematic Biology* 59 (3), 307-321, 2010.
- Accounting for calibration uncertainty: Bayesian molecular dating as a “doubly intractable” problem S Guindon *Systematic Biology* 67 (4), 651–661, 2018.
- Demographic inference under the coalescent in a spatial continuum S Guindon, H Guo, D Welch *Theoretical Population Biology* 111, 43-50, 2016.


### Installation

To install any program that is part of the PhyML package, type the following command:

```bash
sh ./autogen.sh;
```

If you are using a Mac computer or running a Unix-like operating system, you will need 
to install the packages autoconf automake and pkg-config. On a Mac, the following command 
should set you up (provided Homebrew is installed on your Mac...): brew install pkg-config
autoconf automake;

Next, to install any program that is part of the PhyML package, type the following commands:

```bash
./configure --enable-XXXX;
make;
```
where XXXX is phyml or phyrex or phytime.

To compile a Windows executable, install MinGW and run:

```bash
./configure --enable-win --enable-XXXX;
make;
```

To install the MPI version of PhyML, type the following commands:

```bash
autoreconf -i;
./configure --enable-phyml-mpi;
make;
```

### Multithreading

`configure` enables OpenMP when the compiler supports it (`--disable-openmp` turns it off), and
`phyml --threads N` then computes the likelihood with N threads. The results are the same for any
number of threads. GCC supports OpenMP out of the box; with Apple's clang, install `libomp`
(`brew install libomp`) and configure with:

```bash
./configure --enable-phyml \
  CPPFLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
  LIBS="-L$(brew --prefix libomp)/lib -lomp"
```

### Tests

`make -C test check` builds its own binaries and runs the C unit tests and the Python functional
tests (a few minutes; only a C compiler and Python 3 are needed). `make -C test check-mpi` tests the
MPI bootstrap. See `test/README.md`. GitHub Actions and GitLab CI run them on every push.

### Documentation

The user manual is `doc/phyml-manual.pdf`; `phyml --help` lists all the options.
