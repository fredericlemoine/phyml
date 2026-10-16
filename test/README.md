# PhyML test suite

Tests for the `phyml` program as it is in this repository. Nothing outside a C
compiler, `make` and Python 3 (standard library only) is needed.

```sh
make -C test check            # build test binaries, run unit + functional tests
make -C test check-unit       # C unit tests only
make -C test check-functional # Python tests only (TESTS='test_cli*' to filter)
make -C test check-asan       # whole functional suite on the sanitizer build (slow)
make -C test check-mpi        # MPI bootstrap tests (needs mpicc and mpirun)
```

The suite builds its own binaries from `src/` into `test/build/`, so it does
not depend on autotools or on a previous `./configure`:

| Binary | Flags | Used for |
| --- | --- | --- |
| `phyml` | `-O2`, plus `-march=native` on x86-64 (as `configure`) | everything |
| `phyml-generic` | `-DDISABLE_NATIVE` | SIMD kernels vs scalar code |
| `phyml-asan` | `-fsanitize=address,undefined -fsanitize-recover` | memory-safety tests |
| `phyml-ndebug` | `-DNDEBUG` | behaviour without assertions |
| `phyml-omp` | OpenMP (`-fopenmp`; Homebrew `libomp` on macOS), built when the compiler supports it | `--threads` |
| `phyml-mpi` | `mpicc -DMPI` (only for `check-mpi`) | MPI bootstrap |

Useful variables: `CC`, `OPTFLAGS`, `ARCHFLAGS`, `PYTHON`, `SANITIZE=0` (skip
the sanitizer build), `OPENMP=0` (skip the OpenMP build), `OPENMP_CFLAGS` and
`OPENMP_LIBS`, `MPICC`, `MPIRUN` (e.g. `MPIRUN="mpirun --oversubscribe"`).
To test another binary, run the Python tests with `PHYML=/path/to/phyml`.

## Layout

| Path | Content |
| --- | --- |
| `unit/` | C unit tests linked against PhyML's objects (`tinytest.h`: each test runs in a forked child, so `exit()` or a crash fails one test only) |
| `functional/` | Python `unittest` tests that run the `phyml` binary |
| `functional/oracle.py` | Independent reference implementation of the likelihood (Felsenstein pruning, GTR family, discrete Gamma, +I, ambiguity codes, custom amino-acid matrices) and of marginal ancestral states |
| `functional/newick.py` | Newick parser, splits, Robinson-Foulds and transfer-bootstrap reference code |
| `data/` | Small alignments and trees; regenerate with `make -C test data` (`tools/make_data.py`, deterministic) |
| `data/aa_matrices/` | PAML files of protein matrices tested by SMS: Q.pfam, Q.bird, Q.insect, Q.mammal, Q.plant, Q.yeast and MtZoa (from SMS 2.1.1, identical to IQ-TREE's), AB (identical to the paper's supplementary file) and FLU; not generated |

## What is tested

| Area | Unit (C) | Functional (Python) |
| --- | --- | --- |
| Random numbers, distributions, discrete Gamma, aLRT conversion | `test_stats.c` | |
| Newick reading and writing, bipartitions, RF, transfer distance, BioNJ | `test_trees.c` | |
| State encoding, strings, matrix inverse | `test_misc.c` | |
| Likelihood vs. the reference implementation (JC69, K80, F81, HKY85, GTR, +G mean and median, +I, +I+G, ambiguities, custom AA matrix, the protein matrices added for SMS, per-site values, SIMD vs scalar) | | `test_likelihood_oracle.py` |
| Models and optimisation (nesting, parameter recovery on simulated data, `-o` levels, fixed values, FreeRate, invariance to taxon and site order) | | `test_models.py` |
| Tree search (true topology recovered, optimality, starting trees, random starts, constraints) | | `test_search.py` |
| Supports (aLRT, chi2, SH, aBayes, bootstrap counts, TBE vs reference, Bayesian bootstrap, reproducibility) | | `test_supports.py` |
| Input formats, duplicates, malformed alignments and trees | | `test_inputs.py` |
| Command line, output names, seeds, multiple data sets | | `test_cli.py` |
| Ancestral states vs reference, site weights, traces, stats file | | `test_outputs.py` |
| Memory errors under ASan/UBSan | | `test_memory_safety.py` |
| `-DNDEBUG` build | | `test_build_variants.py` |
| Multithreading: identical results with 1, 2 and 4 threads (search, supports, bootstrap, ancestral states, site likelihoods), `--threads` validation | | `test_threads.py` |
| SMS model selection (`--sms`): fits table (K, AIC, BIC, order), heuristic path, nested models, selected model used by the analysis, `--sms_only`, `--sms_criterion`, `-c`, `-o tl`, `-u`, several data sets, refused options, same fits for any number of threads | | `test_sms.py` |
| MPI bootstrap | | `test_mpi.py` |

## Known bugs

Tests that check behaviour PhyML does not have yet are marked as known bugs
(`known_bug(...)` in Python, `KNOWN_BUG(...)` in C). Each names the defect and,
where known, the source line. They are reported as expected failures
(`x` / `XFAIL`), so `make check` passes on the current code. When a bug gets
fixed, its test reports an unexpected success (`XPASS`) and the run fails until
the marker is removed.

To see how each known-bug test currently fails, run the Python tests with
`PHYML_SHOW_KNOWN_BUGS=1`, or run `build/unit_tests` and read the `XFAIL`
lines.

Known bugs covered at the time of writing:

- Crashes: `-d generic` aborts whenever parameters are optimised
  (`eigen.c:53`); the covarion model (`--cov`) segfaults; internal
  polytomies (`io.c:247`); unbalanced or
  unclosed Newick input; truncated or short sequences and fewer than three taxa
  end in `assert()`.
- Memory errors found by ASan: `%1000s` taxon-name reads (`Read_Seq_Sequential`/`Read_Seq_Interleaved`),
  1000-character names, sequence lines longer than the header (`io.c:1721`),
  long input paths (`cl.c:1842`) and the "file name too long" message itself
  (`cl.c:1415`).
- Wrong results: constraints not respected by the search; root polytomies
  resolved with copied branch lengths; with `-n`, random starts only on the
  first data set and two trees per data set in the tree file; MPI silently
  rounds the replicate count up;
  `Permutate()` is not a uniform shuffle.
- Other: `-DNDEBUG` builds hang on sequential files (reader inside `assert`);
  `--print_trace` writes an empty file; duplicate taxon names are accepted;
  `--version`/`--help` exit with status 1; stale `config.h` version; help text
  out of date (`-s NNI`, `-b` range).

## Notes

- PhyML's `-t` is used as κ (the transition/transversion *rate* ratio) for K80
  and HKY85, not as the ts/tv ratio described in `--help`; the oracle tests
  encode the actual behaviour.
- Searches and bootstraps use `rand()`, whose sequence differs between C
  libraries, so search tests check properties (true topology, optimality,
  ranges, invariants) rather than exact numbers. Exact values are only
  compared for fixed-parameter likelihoods, and between runs of the same
  binary (thread counts, reruns) with a fixed `--r_seed`: without it the seed
  comes from the clock.
- On x86-64, `src/sse.c` does not compile without `-msse3` or
  `-march=native`, because the SSE kernels are enabled by `__SSE2__` but use
  SSE3 intrinsics. The test Makefile uses `-march=native` there, like
  `configure`.
