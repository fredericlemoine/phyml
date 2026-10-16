# PhyML: Software Report and Code Assessment

2026-10-14 · Frederic Lemoine (section 3 updated with test-suite results; multithreading and SMS model selection added 2026-10-15)

## 1. Introduction

PhyML is a C program that infers phylogenetic trees from aligned DNA or protein sequences under the maximum-likelihood (ML) criterion. The repository (version 3.3.20260528, about 132,000 lines of C in 47 modules, about 2,100 commits since 2014, almost all by Stéphane Guindon) also ships three Bayesian companion tools built from the same code base.

**Domain.** Molecular phylogenetics reconstructs the evolutionary history of genes, species or pathogens from sequence data. ML methods score a tree by the probability of the alignment given the topology, branch lengths and a Markov model of substitution, then search for the tree that maximises it. The number of unrooted topologies grows super-exponentially (about 2 × 10^182 for 100 taxa), so practical tools pair fast likelihood evaluation with heuristic topology moves.

**Goals.** PhyML was published in 2003 (Guindon & Gascuel) as "simple, fast and accurate": reach likelihoods close to those of slow exhaustive searches by optimising topology and branch lengths together. Version 3.0 (Guindon et al., Syst. Biol. 2010) added SPR moves and the aLRT branch tests. It is among the most cited tools in biology and powers the ATGC Montpellier web server.

**Programs in the package** (one binary per `./configure --enable-XXX`, selected by preprocessor macros):

| Binary | Purpose | Framework |
| --- | --- | --- |
| `phyml`, `phyml-mpi` | ML tree, model parameters, branch supports (bootstrap, TBE, aLRT, aBayes) | Maximum likelihood |
| `phytime` | Divergence dating with calibration uncertainty and relaxed clocks | Bayesian MCMC |
| `phyrex`, `phyrexsim` | Spatial Lambda-Fleming-Viot coalescent: population density and dispersal from geo-referenced sequences | Bayesian MCMC, simulation |
| `evolve`, `date`, `rf`, `test` | Sequence simulation, dating, Robinson-Foulds distances, internal tests | Utilities |

The audience is biologists running single analyses from the command line, an interactive PHYLIP-style menu or an XML file, and pipelines that call `phyml` in batch (it is packaged in Bioconda and Debian).

## 2. Algorithms and comparison with competitors

PhyML combines a classical pruning-based likelihood engine with a staged, perturbed SPR search and offers the widest menu of fast branch-support tests in the field; its weak spots are scalability (high memory, and multithreading only just added) and search depth on large data sets.

### 2.1 Likelihood engine (`lk.c`, `sse.c`, `avx.c`, `eigen.c`)

- **Felsenstein pruning** with partial likelihoods stored on both sides of every edge (`p_lk_left`, `p_lk_rght`), filled by a post-order then a pre-order traversal. Any edge can then be scored in O(sites) without a new traversal, which makes branch-length optimisation and SPR scoring cheap.
- **Transition matrices** from an eigen-decomposition of the rate matrix (`eigen.c`), one per edge and rate class; an "eigen-LR" variant precomputes dot products so a branch length change costs only an `exp()` per state.
- **Numerical scaling**: integer scaling exponents per site (`sum_scale`), corrected in `Rate_Correction`.
- **SIMD**: hand-written SSE and AVX kernels for 4 states (DNA) and 20 states (protein), selected at compile time; ARM goes through the `sse2neon.h` translation layer. Other state counts use a scalar generic path.
- **Site-pattern compression** (`Compact_Data`) and optional sub-pattern aliasing (`--alias_subpatt`). An optional BEAGLE back end exists but rejects bootstrap and +I models.

### 2.2 Models

- **DNA**: JC69, K80, F81, F84, HKY85, TN93, GTR, any custom 6-letter rate grouping (`-m 012345`) and user-given GTR rates (both crashed until 2026-10-16, section 3.1). `-t` sets kappa, the transition/transversion rate ratio.
- **Protein**: LG, WAG, JTT, Dayhoff, DCMut, CpREV, MtREV, MtMam, MtArt, RtREV, VT, Blosum62, HIVb, HIVw, AB, FLU, plus user matrices and the LG4M and LG4X mixtures. Q.pfam, Q.bird, Q.insect, Q.mammal, Q.plant, Q.yeast (QMaker) and MtZoa were added with SMS (below); their likelihoods match the reference implementation computed from the PAML files.
- **Rate heterogeneity**: discrete Gamma (+G), invariant sites (+I), FreeRate (+R), the covarion-like M4 model (`m4.c`), and mixtures or partitions described in XML (`mixt.c`, `xml.c`).
- **Model selection (new)**: `--sms` runs SMS (Smart Model Selection, Lefort et al. 2017) inside PhyML before the analysis, which then uses the selected model; `--sms_criterion aic|bic` and `--sms_only` complete it (`sms.c`). SMS was a set of bash scripts calling the PhyML binary once per candidate, plus R and `bc`. The heuristic is unchanged: fixed BioNJ (or `-u`) tree with GTR+G+I or LG, then the best rate decoration, the best matrix (GTR > TN93 > HKY85 > K80 for DNA; the 23 protein matrices ordered by chi-square distance of their frequencies, with and without +F), and the final decoration. On four test alignments (DNA; small and large protein; protein with a p-inv fixed at Step 3) it gives the same fits as the SMS 2.1.1 scripts run with the same PhyML, within 1e-5 in log-likelihood, and the same selected model on three of them. On the fourth (`aa_small`), the scripts select Q.yeast+G+I although Q.yeast+G has the better BIC: their Step 4 for proteins (`bestFinalAA` in `strategy.sh`) computes the +G fit but never compares it, so +G alone is never selected for proteins. `--sms` compares all the decorations. A second difference: p-inv is passed at full precision rather than rounded to 3 decimals (up to 7e-3 lnL). With the unmodified SMS 2.1.1 head (commit 916269c, run on Linux with R 4.2.3, since its `chi.R` fails with R 4.3 and later), the results are the same, up to the noise of its unseeded runs; with its bundled PhyML 3.3.20220809, FreeRate fits differ by up to 5.5 lnL, but the selected models are unchanged. The built-in FLU matrix now has the full precision of SMS's file (it was rounded to 6 decimals, a 3e-4 lnL difference). Fits are written to `*_phyml_sms.csv`.
- No codon models, no profile mixtures such as C60 or PMSF.

### 2.3 Parameter optimisation (`optimiz.c`)

Brent line searches and a Newton-Raphson/spline method for branch lengths (`Br_Len_Opt`, `Br_Len_Spline`), BFGS for groups of parameters, coordinate-wise rounds over all free parameters (`Round_Optimize`). Pairwise ML distances feed BioNJ (`bionj.c`) for the starting tree.

### 2.4 Tree search (`spr.c`, `Global_Spr_Search`)

The default search is the same whatever `-s` says (NNI and BEST are deprecated and silently mapped to it):

1. Build a BioNJ tree from ML distances (or read a user, random or constraint tree) and optimise model parameters.
2. **Phase 1**: repeated rounds of SPR regrafting with a large radius (max(40, n/4) edges), each preceded by a random Gamma perturbation of all branch lengths whose variance shrinks over time, an annealing-like escape from local optima. Moves are filtered by a quick likelihood bound (`max_delta_lnL_spr`) and applied greedily.
3. **Phase 2**: tighter radius (at most 20) and tolerances until fewer than 5 improvements per round.
4. **Phase 3**: radius at most 30, candidate regraft lists re-evaluated with "triple" moves (local branch-length re-optimisation around the move), stop after 10 rounds without a better tree.
5. Restore the best tree seen, then alternate full parameter rounds and five-branch NNI checks until no NNI improves.

Parsimony-guided SPR (`Spr_Pars`), random starting trees and multiple random starts are available as options.

### 2.5 Branch supports

| Method | Option | Cost | Notes |
| --- | --- | --- | --- |
| Non-parametric bootstrap (Felsenstein) | `-b N` | N full searches | MPI parallel over replicates |
| Transfer Bootstrap Expectation | `--tbe -b N` | N full searches | Lemoine et al., Nature 2018; better for deep branches of large trees |
| Bayesian bootstrap | `--bayesian_bootstrap -b N` | N full searches | Dirichlet site weights |
| aLRT statistic, chi2-based aLRT | `-b -1`, `-b -2` | One NNI pass | Anisimova & Gascuel 2006 |
| min(chi2, SH-like) | `-b -3` | One NNI pass + RELL | Conservative combination |
| SH-like aLRT | `-b -4` | One NNI pass + 10,000 RELL resamples per branch | Guindon et al. 2010 |
| aBayes | `-b -5` (default) | One NNI pass | Anisimova et al. 2011 |

### 2.6 Bayesian tools

- **PhyTime** (`times.c`, `rates.c`, `mcmc.c`): MCMC dating under autocorrelated and uncorrelated relaxed clocks, with calibrations treated as uncertain clade placements, a "doubly intractable" posterior (Guindon 2018). A normal approximation of the likelihood around the ML branch lengths makes it fast.
- **PhyREX** (`phyrex.c`, `slfv.c`, `rrw.c`, `ibm.c`, `velocity.c`): spatial Lambda-Fleming-Viot coalescent, relaxed random walks and integrated velocity models of dispersal, output readable by TreeAnnotator and SPREAD.

### 2.7 Comparison with competitors

| Tool | Search strategy | Speed and scale | Models | Supports | Parallelism |
| --- | --- | --- | --- | --- | --- |
| **PhyML 3.3** | BioNJ + perturbed SPR + NNI polish | Fine up to a few thousand taxa; hard limit of 4,000 taxa | Rich DNA/AA, +G +I +R, LG4X, XML mixtures; no codons, no model test | Bootstrap, TBE, Bayesian bootstrap, aLRT, aBayes, SH-aLRT | MPI for bootstrap; OpenMP over sites (`--threads`, new) |
| **RAxML-NG** | Lazy SPR with adaptive radius | Very fast (libpll, AVX2), 10^4 to 10^5 taxa | DNA/AA, binary, morphology, partitions, +G +I +R | Bootstrap, TBE, convergence test | Pthreads + MPI, fine-grained over sites |
| **IQ-TREE 2/3** | Stochastic NNI with a population of candidate trees | Fast, usually highest likelihoods in benchmarks | Largest model set incl. codons, profile mixtures, ModelFinder | UFBoot, SH-aLRT, aBayes, concordance factors | OpenMP + MPI |
| **FastTree 2** | Minimum evolution + NNI/SPR + ML NNI | Fastest; 10^5 to 10^6 taxa | Limited (JTT, WAG, LG, GTR, CAT approximation) | SH-like local supports | OpenMP |
| **MrBayes, BEAST 2** (vs PhyTime, PhyREX) | MCMC | Slower; full posterior | Large Bayesian model space, many clocks and phylogeography models | Posterior probabilities | BEAGLE, MC3 |

Independent benchmarks on empirical phylogenomic data (Zhou et al., Mol. Biol. Evol. 2018) found IQ-TREE and RAxML reaching higher likelihoods than PhyML in most cases, with PhyML often slower on large alignments. PhyML's distinctive strengths are the fast local supports it pioneered (aLRT, aBayes, SH-aLRT, now copied by IQ-TREE), TBE, its well-validated defaults for single-gene analyses, and a long record of stability. The comparison above comes from published literature; this review did not rerun benchmarks.

## 3. Critical assessment of the code

The numerical core is sound: on fixed trees, PhyML's likelihoods match an independent reference implementation to 1e-8 for every model tested. The code around it carries real defects, though: 29 of the 35 problems listed below were reproduced by running PhyML, 11 of them crashes, aborts or memory errors. Four are now fixed: `-o n` without `-m`, custom DNA models, the AB matrix, and TBE supports in `phyml-mpi`. On top of these come unchecked buffers, no automated tests until now, and a monolithic structure that makes each fix risky. The first version of this report came from reading the source at commit 42bb4b5. The findings were then checked with the test suite in `test/` (section 3.4), run on macOS arm64, Linux aarch64 and Linux x86-64, with sanitizer, `-DNDEBUG` and MPI builds.

### 3.1 Bugs and correctness risks

"Verified" says how each item was checked: **test** = reproduced by a known-bug test in `test/` (file named), **ASan** = reproduced by the AddressSanitizer build, **code** = found by reading only.

| Severity | Location | Problem | Consequence | Verified |
| --- | --- | --- | --- | --- |
| High | `mpi_boot.c` (~line 330), `io.c:6818` | `phyml-mpi` reduces only `bip_score` across ranks. TBE's `tdist_score` stays local to rank 0 but is divided by the total replicate count. | `phyml-mpi --tbe` reported supports biased towards 1: on a test branch it printed 1.0 where its own replicate trees give 0.9. | Fixed 2026-10-16 (`tdist_score` is now reduced too); test (`test_mpi.py`) |
| High | `cl.c:956-968`, `cl.c:980` | Custom DNA models (`-m 012345`, `-m 000000`, ...) and user GTR rates (`-m 1.5,3,0.7,1.2,4,1`) write to `io->mod->r_mat` before it is allocated. | Segmentation fault for every custom nucleotide model, although `--help` documents them. | Fixed 2026-10-16, with three related defects: user GTR rates were then re-estimated (reset by `Init_Model`), custom models with `-f o` started from all-zero frequencies, and the stats file printed garbage as the custom model name; test (`test_likelihood_oracle.py`, `test_models.py`) |
| High | `cl.c:1882` | With `-o n` and no `-m`, the option handler writes `io->mod->r_mat->optimize` while `r_mat` is still NULL. | `phyml -i aln -o n` crashes (also `-o tl`, `-o l`); it works only when `-m` is given. | Fixed 2026-10-15 (needed by `--sms`, which forbids `-m`); test (`test_cli.py`) |
| High | `io.c:327` (`Read_Node_Name`) | Copies a tip name until `:` or `[`; a tip without a branch length has neither. | Heap overflow and crash on any user tree without branch lengths, and therefore on every `--constraint_file` (constraint trees have no lengths). | test, ASan (`test_inputs.py`, `test_search.py`, `test_memory_safety.py`) |
| High | `io.c:247-256` (`R_rtree`) | After `Clean_Multifurcation`, frees the subtree string, allocates an empty buffer and parses it again without copying the resolved subtrees into it. | Any user tree with an internal polytomy crashes. | test (`test_inputs.py`, `test_trees.c`) |
| High | `eigen.c:53`, generic data type | With `-d generic`, the eigen-decomposition of the rate matrix produces NaN as soon as parameters are optimised. | Every generic-data analysis except `-o n` aborts with "Assertion failed", for 2, 3 or 5 states alike. | test (`test_models.py`) |
| High | Covarion model (`--cov`, `m4.c`) | The run segfaults right after the first optimisation of the switching rate lambda. | The covarion model, documented in `--help`, cannot be used. | test (`test_models.py`) |
| High | `main.c:355` | Every rank runs its own full ML search (random Gamma perturbations, seed + rank), then `MPI_Bcast` uses each rank's own `strlen` as count. | Mismatched counts are undefined in MPI; a receiver whose tree string is shorter than rank 0's can be overrun. All ranks also redo the whole ML search for nothing. | code |
| High | `stats.c:88` | `Rand_Int` returns `rand() % n`. On the Windows build `RAND_MAX` is 32,767. | Windows bootstrap replicates on alignments longer than 32,768 sites never sample the remaining sites; modulo bias elsewhere. | code (a unit test passes on Unix, where `RAND_MAX` is 2^31-1) |
| High | `io.c:1675-1725` | `Read_One_Line_Seq` writes `state[len++]` with no bound, and stores `fgetc` in a `char` before comparing with `EOF`. | Heap overflow when a line holds more characters than the header declares (ASan: `io.c:1721`). Where `char` is unsigned (Linux aarch64) EOF is never seen. | ASan (`test_memory_safety.py`) |
| Medium | `init.c` (`Init_Qmat_AB`) | The built-in AB matrix swaps the Thr-Arg and Thr-Asn exchangeabilities (1.135 and 2.147) of the `AB_model` file published with the paper (Mirsky et al. 2015, supplementary material 1). | `-m AB` gave wrong likelihoods (0.3 to 0.8 lnL units on the test alignments). | Fixed 2026-10-15; test (`test_likelihood_oracle.py`, against the published file) |
| Medium | Constraint search (`spr.c`, `main.c:180-205`) | With a constraint naming all taxa (the only form that does not crash), the search moves away from it. | `((t1,t3),...)` as constraint, output tree groups t1 with t2: the constraint is silently ignored. | test (`test_search.py`) |
| Medium | `io.c:1504` | `assert(Read_One_Line_Seq(...))` puts the read inside an assert. | A `-DNDEBUG` build hangs forever on sequential files; today short or truncated files abort with an assert instead of a message. | test (`test_build_variants.py`, `test_inputs.py`) |
| Medium | `io.c:1487`, `io.c:1548` | Format `%1000s` reads into a 1,000-byte name buffer. | One-byte overflow for 1,000-character taxon names (ASan flags it on every run); `Check_Sequence_Name` overflows too (`io.c:1501`). | ASan (`test_memory_safety.py`) |
| Medium | `cl.c:267-285` and four similar blocks, `cl.c:1842`, `interface.c:128-192` | The "file name too long" error path `strcat`s the over-long name into a 200-byte buffer; output names append suffixes plus a run id to 200-byte buffers that accept 189-character paths. | Heap overflow on long paths, common on HPC scratch directories (ASan: `cl.c:1415`, `cl.c:1842`). | ASan (`test_memory_safety.py`) |
| Medium | Newick reader (`io.c:39`, `Sub_Trees`) | An unclosed `[` makes `Read_Tree` scan past the string; unbalanced parentheses are not detected; `[comments]` are skipped when counting taxa but not when parsing. | Bus error on malformed trees instead of an error message; valid Newick with comments is rejected. | test, ASan (`test_inputs.py`, `test_trees.c`) |
| Medium | `Clean_Multifurcation` (`io.c`) | Edges created to resolve a root polytomy get a copy of a sibling's length instead of zero. | `(A:.1,B:.2,C:.3,D:.4,E:.5)` is read with total length 2.2 instead of 1.5. | test (`test_trees.c`) |
| Medium | `main.c:331-334` | After the first data set, `random_input_tree` is set to NO and never restored. | With `-n` > 1 and `--rand_start`, data sets 2 to n silently start from BioNJ only. | test (`test_cli.py`) |
| Medium | `main.c:384-397` | With `-n` > 1, the tree is written before and after re-inserting duplicates, without rewinding. | `*_phyml_tree.txt` holds two trees per data set (6 for 3 data sets). | test (`test_cli.py`) |
| Medium | `utilities.c:8166-8250` | `Check_Memory_Amount` multiplies `int`s before widening; asks an interactive question above 256 MB. | Overflowed (even negative) estimates for large data; batch jobs wait on stdin unless `--no_memory_check` is set. | code |
| Medium | `mpi_boot.c:76-82`, `:255-300` | Replicate count silently rounded up, message printed by every rank; each tree sent as a fixed 2 MB `MPI_Ssend`, received rank by rank. | 10 replicates on 3 ranks gave 12; communication cost unrelated to tree size; trees longer than 2 MB are truncated. | test (`test_mpi.py`) |
| Medium | `io.c:1105`, `io.c:1504` | Fewer than three taxa, short sequences and truncated files end in `assert()`. | Abort with "Assertion failed" instead of an error message. | test (`test_inputs.py`) |
| Medium | Sequence reader | Duplicate taxon names are accepted. | The analysis runs and the output tree has two tips with the same name. | test (`test_inputs.py`) |
| Low | `help.c:245`, `cl.c:862-882` | Help says NNI is the fast default, but NNI and BEST map to the same SPR search, and `-s NNI` silently turns off random starts. | Users and wrappers rely on options that do nothing. | test (`test_cli.py`, `test_search.py`) |
| Low | `help.c` | `-t` is documented as the ts/tv ratio but is used as kappa (the rate ratio): K80 with `-t 2` matches kappa = 2 exactly. | Users converting between programs set the wrong value. | test (`test_likelihood_oracle.py`) |
| Low | `--print_trace` | The option creates `*_phyml_trace.txt` but writes nothing to it. | Documented feature missing (the JSON trace works). | test (`test_outputs.py`) |
| Low | `utilities.c:1115` | `mCalloc` takes an `int` count, treats `calloc(0)` returning NULL as failure, and handles out-of-memory with `assert(false)`. | With `-DNDEBUG`, allocation failure returns NULL and crashes later. | code |
| Low | `utilities.c:3990-4010` | Bootstrap with a user tree rewinds the tree file and reads its first tree; replicates switch off `print_site_lnl` on the shared `io`. | Wrong starting tree with multi-tree files; hidden side effects on later output. | code |
| Low | `cl.c:1198-1245` | `-b` error message says -1 to -4 while -5 (aBayes, the default) is valid; the bootstrap versus multiple data sets check depends on option order. | Confusing errors; `-b 100 -n 5` and `-n 5 -b 100` behave differently. | test (message only, `test_cli.py`) |
| Low | `stats.c:4659` | `Permutate` swaps `x[i]` with a random position over the whole array instead of Fisher-Yates. | Biased permutations (for 3 elements, some orders have probability 4/27, others 5/27); used to order SPR candidates. | test (`test_stats.c`) |
| Low | `stats.c:76` | `Uni()` returns values in [0, 1] inclusive. | Possible `log(0)` in samplers that use `log(u)` or `log(1-u)` (PhyTime, PhyREX). | code |
| Low | `main.c`, `help.c` | `--version` and `--help` exit with status 1. | Scripts that probe the version see a failure. | test (`test_cli.py`) |
| Low | `config.h` | The committed `config.h` says 3.3.20250429 while `configure.ac` says 3.3.20260528. | Builds that do not rerun `configure` report the wrong version. | test (`test_cli.py`) |
| Low | `utilities.h:37-45`, `sse.c` | On x86-64 the SSE kernels are enabled by `__SSE2__` but call the SSE3 intrinsic `_mm_hadd_pd`. | PhyML does not compile with the compiler's default flags; `configure` hides it by always adding `-march=native`. | build (Linux x86-64) |
| Low | `Dist_Btw_Tips` (`utilities.c:11649`) | `Path_Length` relies on `t_dir` routing tables that `Read_Tree` does not set up. | Calling it on a freshly parsed tree loops forever (found while writing unit tests). | unit-test probe |

Behaviour that was checked and found correct includes all fixed-parameter likelihoods, per-site likelihoods, ancestral marginal probabilities, TBE in serial runs, bootstrap counts, the effect of `--weights`, identical results from the SIMD and scalar kernels on x86-64, and identical results for any number of threads.

The history shows the cost of fixes without tests: 361 of about 2,100 commits mention a fix, several repeat the same subject ("Fixed conceptual bug in kfol-pos method" four times), and commit dc44563, "Sanity checks on alrt routines", changes only generated build files.

### 3.2 Performance

- **Shared-memory parallelism: first step done.** Until this update the only parallelism was MPI across bootstrap replicates. `--threads N` now splits the site patterns of every likelihood, derivative, partial-likelihood and eigen-projection loop over OpenMP threads (`lk.c`, `sse.c`, `avx.c`). Each thread uses its own scratch vectors, and per-site values are summed in site order, so results are bit-identical for any thread count. On a 40-taxon, 20,000-site GTR+G search on an Apple M-series laptop (4 performance cores), the run took 59.9 s with 1 thread, 35.1 s with 2 and 22.9 s with 4 (2.6x). Remaining serial parts: mixture and partition models (`mixt.c`), parsimony, distance matrices and per-branch support tests. Threads only engage above 128 patterns per thread, so short alignments gain little. The default stays at 1 thread.
- **Memory layout.** Conditional likelihoods are stored on both sides of each of the 2n-3 edges, about twice a node-based layout. Parsimony vectors use a 4-byte `int` per state per site instead of bit vectors. This, plus the hard cap `N_MAX_OTU` = 4,000 taxa, limits PhyML to mid-sized data.
- **Likelihood loop.** Every full evaluation redoes a pre-order pass even when only one edge is scored next.
- **SIMD at compile time only.** `configure.ac` defaults to `-march=native` and picks SSE or AVX by macro. Distributed binaries either crash with an illegal instruction on older CPUs or fall back to the scalar path; there is no runtime dispatch and no AVX-512 kernel.
- **SH-like aLRT.** `Statistics_To_SH` allocates and frees `pi` and a sample array inside its 10,000-iteration loop, for every branch, although `pi` is constant.
- **Bootstrap replicates** recompute an ML distance matrix and BioNJ tree (O(n^2 L)), deep-copy the model, and for TBE allocate three n x n matrices per replicate; `Match_Tip_Numbers` is an O(n^2) `strcmp` loop per replicate. Faster TBE algorithms exist (Truszkowski, Gascuel & Swenson 2020).
- **Search heuristics.** `Global_Spr_Search` re-scales every branch length by a Gamma draw before each SPR round; the `freq` counter meant to space these out is always 1 (dead logic), so each round starts by undoing part of the previous optimisation.

### 3.3 Architecture and maintainability

- **Size and coupling.** 132,000 lines; `utilities.c` (14,877 lines) and `mcmc.c` (14,901) are catch-all files, and `utilities.h` (2,986 lines) defines every struct. `t_tree` carries ML, MCMC, dating and spatial fields at once; every binary compiles most modules.
- **Very long functions.** `Read_Command_Line` is 1,916 lines, `XML_Process_Base` 1,220, `EVOLVE_Main` 1,049, `Launch_Interface_Model` 882.
- **Dead and rotting code.** About 4,200 lines of commented-out code, 462 preprocessor conditionals in `.c` files, and four modules outside the build (`geo.c`, which no longer compiles, `invitee.c`, `checkpoint.c`, `beagle_utils.c`). Since `beagle_utils.c` is not in `Makefile.am`, the BEAGLE back end cannot be built from the distributed scripts.
- **Aliasing by design.** Bootstrap trees share `s_opt`, `io` and likelihood buffers with the main tree (`Share_Lk_Struct`, comments marked WARNING), so any write in one context leaks into another.
- **Error handling.** About 1,000 calls to `exit`, `Exit` or `assert(FALSE)` deep in the code, and 34 "Type any key to exit" prompts that wait on stdin. The engine cannot be embedded as a library (the R wrapper is commented out).
- **Type hygiene.** 97 `-Wsign-compare` and 121 unused-parameter warnings; `unsigned int` variables set to -1 in `Lk()`; only 4 `static` functions in the whole code base.

### 3.4 Build, packaging and testing

- **Tests: now a first suite, not yet in CI.** The repository had no tests, and GitHub Actions builds on three platforms without running an analysis. `test/` now holds 56 C unit tests and 196 Python functional tests (`make -C test check`, about 1 to 2.5 minutes). They compare likelihoods with an independent Python implementation, check search, support and optimisation properties on simulated data, and run sanitizer, `-DNDEBUG`, OpenMP and MPI builds. The 34 known-bug tests (33 in `make check`, 1 more in `make check-mpi`) are reported as expected failures, so the suite passes on the current code and flags each fix. The README still shows a Travis badge (travis-ci.org closed in 2021).
- **Generated files in git.** `config.h`, `config.log`, `config.status`, the generated `Makefile`s and LaTeX `.aux`/`.log` files are tracked; each commit carries about 1,000 lines of build noise, which hides real changes in review.
- **Configure quirks.** `AC_DEFINE([DEBUG],[1])` is unconditional; the compiler-specific `-mcpu=native` choice is overwritten a few lines later; one binary per configure run, chosen by ten mutually exclusive flags. The unconditional `-march=native` also hides that the x86-64 SSE code needs SSE3 (section 3.1).
- **Sanitizer cost.** Under GCC 13's AddressSanitizer, PhyML runs about 40 times slower unless `detect_stack_use_after_return=0` is set; the test suite sets it.

## 4. Axes for improvement

Build a regression safety net first, then fix the correctness bugs, then invest in parallelism and modularity; the order matters because the later steps are large refactors that only tests make safe.

1. **Regression tests and CI (started).**
    - Done: `test/` holds unit, functional, sanitizer, `-DNDEBUG` and MPI tests, run by `make -C test check` (section 3.4). It passes on macOS arm64, Linux aarch64 and Linux x86-64.
    - Next: run `make -C test check` (and `check-mpi`) in GitHub Actions on all three platforms, so each fix turns its known-bug test into a passing one.
    - Add a fuzzing target for the PHYLIP, NEXUS, Newick and XML parsers; most of the crashes in section 3.1 are in these readers.
    - Remove generated files (`config.*`, `Makefile`s, LaTeX outputs) from git and drop the Travis badge.
2. **Fix the high-severity bugs (days).**
    - Crashes in documented features: make `Read_Node_Name` stop at `,` and `)`; fill the buffer that `R_rtree` re-parses after resolving a polytomy; then check that constraints are kept during SPR moves.
    - MPI: run the ML search on rank 0 only and broadcast the length before the string; keep the user's replicate count and give the remainder to the first ranks.
    - Replace `rand()` with a seeded 64-bit generator (PCG or xoshiro), one stream per MPI rank, unbiased bounded integers and a `Uni()` on the open interval (0, 1).
    - Rewrite the sequence readers with bounds checks, `int` for `fgetc`, real error messages, and no side effects inside `assert`.
    - Replace fixed `T_MAX_*` buffers and `strcpy`/`strcat`/`sprintf` with sized allocations (`snprintf`, `asprintf`); widen size arithmetic to `size_t`.
    - Make PhyML non-interactive by default when stdin is not a terminal; restore random-start state per data set; align the help text with the actual search.
3. **Performance and scale (months).**
    - Threads (started): the site loops of the likelihood engine run on OpenMP threads (`--threads`, section 3.2). Next: mixture and partition models, parsimony and distance matrices, aLRT/aBayes branch tests in parallel, and bootstrap replicates on threads once replicates stop sharing `io` and `s_opt` with the main tree.
    - Replace compile-time `-march=native` with runtime CPU dispatch (scalar, SSE, AVX2, AVX-512, NEON) so one Bioconda binary is both portable and fast.
    - Either adopt libpll-2 or BEAGLE as the kernel (and add `beagle_utils.c` back to the build), or move to node-based conditional likelihoods with bit-vector parsimony to halve memory, then lift the 4,000-taxon cap.
    - Hoist allocations out of hot loops (SH-like resampling, TBE matrices, per-replicate distance matrices) and use a faster transfer-index algorithm for TBE.
    - Re-tune `Global_Spr_Search` against IQ-TREE and RAxML-NG on a public benchmark (for example the Zhou et al. 2018 data sets): fix the dead `freq` schedule, consider keeping a small population of candidate trees.
4. **Architecture (ongoing).**
    - Split the code into a reusable library (`libphyml`: data, models, likelihood, search, supports) and thin front ends for `phyml`, `phytime` and `phyrex`; move `t_tree`'s MCMC and spatial fields into separate structs.
    - Return error codes instead of calling `exit`, so the engine can be embedded (R, Python) and tested function by function.
    - Delete dead modules and commented-out code (git keeps history), break up the 1,000+ line functions, and turn on `-Wextra -Werror` for new code.
    - Replace autotools mode flags with CMake or Meson targets that build all binaries in one pass.
5. **Scientific features to stay competitive.**
    - Model selection (started): SMS is built in (`--sms`). Next: run the independent SMS fits on threads, offer AICc and the FreeRate category count, and add codon models.
    - Profile mixture models (C20 to C60, PMSF) and partition-wise models from a simple NEXUS partition file, not only XML.
    - Checkpoint and restart for long searches and bootstraps (`checkpoint.c` is a start), and a ultrafast-style bootstrap or parametric alternative to make TBE affordable on large trees.

The quickest wins with the highest impact are the Newick reader, the random number generator and the sequence readers. Each is a few dozen lines, and together they remove every high-severity item in section 3.1. The known-bug tests in `test/` show when each one is fixed.
