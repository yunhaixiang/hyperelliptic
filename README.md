# Hyperelliptic Curve Finder

Find presentations of hyperelliptic curves

```text
y^2 = f(x)
```

over odd prime fields `F_p` whose L-polynomial is a trinomial

```text
L(t) = 1 + a_g t^g + p^g t^{2g}
```

allowing the middle coefficient `a_g` to be zero.

The program lists presentations according to the selected reduction mode, not
certified isomorphism classes.

## Usage

```bash
python3 hyperelliptic_finder.py p g [--max N] [--output FILE] [--reduction {pgl2,affine,pgl2save,affinesave}] [--monic-only] [--no-hasse-witt-prefilter] [--quiet|--verbose] [--log]
```

Arguments:

- `p`: odd prime characteristic. Characteristic `2` is not supported for models `y^2 = f(x)`.
- `g`: genus.
- `--max N`: maximum number of presentations to save. Use `0` for the complete search.
- `--output FILE`: output path. Defaults to `hyperelliptic_results.txt`.
- `--reduction {pgl2,affine,pgl2save,affinesave}`: presentation reduction before point counting. Defaults to `pgl2save`.
- `--monic-only`: enumerate only monic presentations. By default, the nonsquare leading-coefficient class is included too.
- Hasse-Witt prefiltering is enabled by default before point counting.
- `--no-hasse-witt-prefilter`: disable the Hasse-Witt prefilter.
- Default output mode: print only accepted presentations, their index, the
  middle coefficient, and presentation/class counts so far.
- `--quiet`: suppress accepted-presentation progress; summaries still print.
- `--verbose`: print every considered presentation and rejection.
- `--log`: also write terminal progress to a matching `.log.txt` file.

If `--output` ends in `.txt`, a matching `.json` file is also written. If it
ends in `.json`, a matching `.txt` file is also written.
By default no progress log file is written.

## Examples

```bash
python3 hyperelliptic_finder.py 3 1 --max 10
python3 hyperelliptic_finder.py 3 2 --output genus2.txt
python3 hyperelliptic_finder.py 5 2 --max 0 --output complete-results.txt
python3 hyperelliptic_finder.py 3 2 --reduction affine --output affine-results.txt
python3 hyperelliptic_finder.py 3 2 --monic-only --output genus2-monic-only.txt
```

## C Version

The repository also includes `hyperelliptic_finder.c`, a standalone C port that
uses FLINT for finite-field irreducible polynomial construction and C
implementations of the search, reduction, Hasse-Witt filtering, and point
counting loops.

Build it with:

```bash
make
```

On this machine the installed Homebrew FLINT library is `x86_64`, so the needed
build command is:

```bash
make ARCH="-arch x86_64"
```

Run it similarly to the Python version:

```bash
./hyperelliptic_finder_c 3 2 --max 10 --output c-results.txt
```

Supported options are `--max`, `--output`, `--reduction
{pgl2,pgl2save,affine,affinesave}`, `--monic-only`,
`--no-hasse-witt-prefilter`, `--quiet`, and `--verbose`. The C binary also
accepts `--workers N --worker-index I`; these are mainly used by `parallel.py`
to run deterministic shards. The C enumerator may find presentations in a
different order than the Python program, so `--max N` does not necessarily
return the same first `N` presentations. Complete runs should agree on totals
and middle-coefficient distributions.

For C batch runs, use `test_c.py`. By default it runs prime `3`, genera `7`
through `20`, and reduction mode `pgl2save`:

```bash
python3 test_c.py --outdir batch_results_c --resume
```

It writes per-case files named like `p3_g7_pgl2save_c.json` and a summary file
`batch_results_c/batch_summary_c.json`. If `hyperelliptic_finder_c` is missing,
the script runs `make ARCH="-arch x86_64"` by default on this machine. Use
`--no-build` if you want it to fail instead. The default per-case timeout is
30 minutes, but a case is only stopped after at least one presentation has been
saved. Use `--case-timeout SECONDS` to change it, or `--case-timeout 0` to
disable it.

## Parallel C Runs

For one C search split across several processes, use `parallel.py`. It runs
`hyperelliptic_finder_c` in deterministic shards and merges the worker JSON
files afterward. By default it uses 6 workers and reduction mode `pgl2save`:

```bash
python3 parallel.py 3 10 --workers 6 --quiet
```

Useful options:

- `--workers N`: number of C worker processes. Defaults to `6`.
- `--outdir DIR`: directory for merged and worker outputs. Defaults to
  `parallel_results`.
- `--output FILE`: merged output path. A matching `.json`/`.txt` file is also
  written.
- `--reduction {pgl2,affine,pgl2save,affinesave}`: reduction mode passed to
  each worker. Defaults to `pgl2save`.
- `--max N`: maximum presentations per worker, not globally. Use `0` for no
  per-worker limit.
- `--dedupe exact`: keep all distinct presentations and remove only exact
  duplicate coefficient lists across workers. This is the default.
- `--dedupe reduction-key`: keep one presentation per reduction class in the
  merged output.
- `--case-timeout SECONDS` or `--timeout SECONDS`: timeout per worker. Defaults
  to `1800` seconds. A worker is only stopped after it has saved at least one
  presentation. Use `0` for no timeout.
- `--resume`: reuse existing worker JSON files.
- `--quiet`, `--verbose`, `--monic-only`, and
  `--no-hasse-witt-prefilter`: forwarded to every worker.

Each worker currently keeps its own matching cache. This means two workers may
independently discover equivalent presentations before the final merge removes
the duplicate data. Parallelism is still useful when point counting dominates,
but using too many workers can waste time on repeated reduction-key and
matching work. A good starting point is 4, 6, or 8 workers, then compare saved
presentations per minute for the target `p` and `g`. The default of 6 is meant
as a practical middle ground, not a universal optimum.

Parallel output includes worker status metadata. If every worker finishes, the
merged file is marked as a complete list. If any worker times out, is missing,
or is interrupted, the merged file is marked partial even though all saved data
from completed or partially completed workers is still merged.

## Batch Runs

Use `test.py` to run many searches and save each case to a separate file.
By default it runs prime `3`, genera `1` through `20`, and reduction mode
`pgl2save`.

```bash
python3 test.py --outdir batch_results --resume
```

Useful options:

- `--outdir DIR`: directory for per-case `.txt`/`.json` output files.
- `--reduction {pgl2,affine,pgl2save,affinesave}`: reduction mode passed to `hyperelliptic_finder.py`. Defaults to `pgl2save`.
- `--max N`: maximum curves per case; `0` means complete search.
- `--case-timeout SECONDS` or `--timeout SECONDS`: maximum runtime per case
  before stopping that run and moving to the next genus/prime, but only after
  at least one presentation has been saved. Defaults to `1800` seconds, i.e. 30
  minutes. Use `0` for no timeout.
- `--min-genus G` and `--max-genus G`: restrict the genus range.
- `--resume`: skip cases whose JSON output already exists.
- Default output mode: stream accepted presentations only.
- `--quiet`: pass `--quiet` to `hyperelliptic_finder.py`.
- `--verbose`: pass `--verbose` to `hyperelliptic_finder.py`.
- `--log`: write `batch.log.txt` and pass `--log` to each finder run.
- `--monic-only`: pass `--monic-only` to each finder run. Output filenames include `_monic`.
- `--no-hasse-witt-prefilter`: disable the default Hasse-Witt prefilter in each finder run.

The script always writes:

```text
batch_results/batch_summary.json
```

When `--log` is passed, the script also writes:

```text
batch_results/batch.log.txt
```

`batch_summary.json` is always updated after each completed case, so partial
progress is preserved if the batch is interrupted with Ctrl+C. Each individual
`hyperelliptic_finder.py` run also saves its own results incrementally. Per-case
progress logs are written only when `--log` is used.

The full grid can be expensive. For a quick test, use:

```bash
python3 test.py --min-genus 1 --max-genus 3 --max 10 --outdir batch_results
```

## LMFDB Comparison

Use `lmfdb.py` to compare local result JSON files with available LMFDB
abelian-variety isogeny-class records over finite fields.

```bash
python3 lmfdb.py p3_g2_results.json --out lmfdb_accuracy.json
```

For batch results, point it at the per-case result files:

```bash
python3 lmfdb.py batch_results/p*_g*.json --out lmfdb_accuracy.json --delay 1.0 --timeout 30
```

Useful options:

- `--out FILE`: output JSON report.
- `--delay SECONDS`: wait between LMFDB requests.
- `--timeout SECONDS`: HTTP timeout per request.

The comparison is by L-polynomial/isogeny class, not by individual curve
presentation. This matters for `pgl2save` and `affinesave`, where many saved
presentations can share the same middle coefficient and therefore the same LMFDB
query. The report includes both `local_presentations` and
`local_unique_l_polynomials`, and for each middle coefficient it lists the saved
presentation indices and canonical presentation indices. LMFDB may not have
records for every genus and prime searched. The reported
`accuracy_among_available` only uses LMFDB records that were found.

If a non-result JSON file such as `batch_summary.json` is passed by accident,
`lmfdb.py` marks it as skipped in the report instead of failing.

## Output

Detailed results are saved to the text and JSON output files. By default, the
terminal prints run status, accepted presentations, their indices, middle
coefficients, and presentation/class counts so far. Use `--verbose` to print
every considered presentation, reduction skip, and rejection. Use `--quiet` to
suppress accepted-presentation progress. By default no progress log file is
written; use `--log` to write the same terminal progress to a `.log.txt` file.

The output metadata includes:

- `monic_enforced`: `false` by default; `true` when `--monic-only` is used.
- `leading_representatives`: `[1, smallest nonsquare]` by default, or `[1]`
  with `--monic-only`.
- `hasse_witt_prefilter`: whether the Cartier-Manin/Hasse-Witt prefilter was used.
- `reduction_class_count`: the number of distinct saved canonical presentation
  indices after the selected reduction.
- `isomorphism_class_count`: the same count for `pgl2` and `pgl2save`; for
  affine modes this is `null` because affine reduction is not full PGL2
  reduction.
- `search_status` and `complete_list`: whether the search finished or was
  interrupted/timed out.
- `total_presentations_found`: number of saved presentations.
- `workers` and `worker_index`: present in C worker output when a sharded run
  is used.

Each saved curve presentation is indexed and includes:

- `f(x)` for the presentation `y^2 = f(x)`.
- The middle L-polynomial coefficient `a_g`.
- The accepted canonical presentation index and polynomial used for reduction.
- In JSON, `f_coeffs`, `canonical_f_coeffs`, and `reduction_key` for scripts
  that need machine-readable matching data.

The text output omits raw coefficient arrays, canonical coefficient arrays, and
machine-only matching flags to keep the result list shorter.

Only the middle coefficient of the L-polynomial is shown because the other
nonzero coefficients in the trinomial form are fixed as `1` and `p^g`.

Progress is saved after every matching presentation. If the program is
interrupted with Ctrl+C, all matches found before the interruption remain in the
output files.

## Reduction Options

The `--reduction` option controls which repeated presentations are skipped
before point counting.

For speed, the search uses a canonical reduction key instead of storing every
orbit member. After a presentation passes the Hasse-Witt prefilter, the program
computes all presentations reached by the chosen reduction mode, normalizes
them, and stores only the lexicographically smallest coefficient tuple as the
dictionary key. Later equivalent presentations have the same key, so they can be
skipped or saved as accepted-equivalent presentations without keeping a large
orbit-member dictionary.

For PGL2 modes, powers of the linear factors used in fractional linear
transformations are precomputed once for each `(p, g)` and reused across orbit
computations. The C version does the same by building a run-level cache of
normalized PGL2 matrices and their linear-factor powers before the search.

The Cartier-Manin/Hasse-Witt prefilter is enabled by default. It computes the
L-polynomial modulo `p` from the Hasse-Witt matrix and rejects a new reduction
orbit before point counting when a pre-middle coefficient is nonzero modulo
`p`. Use `--no-hasse-witt-prefilter` to disable it for comparison or debugging.

### `--reduction pgl2`

This uses fractional linear transformations

```text
x -> (ax + b)/(cx + d)
```

with `ad - bc != 0` over `F_p`. This includes affine transformations and also
transformations such as inversion `x -> 1/x`. It can identify more equivalent
presentations, including some cases where a degree `2g + 1` presentation and a
degree `2g + 2` presentation are related.

Equivalent presentations of an already accepted curve are skipped.

### `--reduction pgl2save`

This is the default. It uses the same PGL2 orbit computation as `pgl2`, but it saves presentations
that are equivalent to already accepted presentations. In that case the program
does not recompute point counts or the L-polynomial; it reuses the known middle
coefficient and records which accepted canonical presentation the saved
presentation is equivalent to.

### `--reduction affine`

This uses only transformations

```text
x -> ax + b
```

with `a != 0` over `F_p`. It is cheaper than `pgl2`, but removes fewer repeated
presentations.

Equivalent presentations of an already accepted curve are skipped.

### `--reduction affinesave`

This uses the same affine orbit computation as `affine`, but it saves
presentations that are equivalent to already accepted presentations. As with
`pgl2save`, it reuses the known middle coefficient and records the accepted
canonical presentation.

By default, transformed equations are normalized to one of two leading
coefficient representatives: `1` for the square class and the smallest
nonsquare in `F_p` for the nonsquare class. With `--monic-only`, transformed
equations are normalized back to monic form only when the corresponding
`y`-scaling exists over `F_p`.

## Notes

The search enumerates squarefree polynomials of degree `2g + 1` and `2g + 2`
over `F_p`. It searches sparse presentations first: leading term only, then
binomials, trinomials, four-term polynomials, and so on by increasing number of
nonzero lower terms. This does not change complete-search coverage, but it often
finds early presentations faster. By default it uses two leading coefficient
representatives: `1` and the smallest nonsquare in `F_p`. This includes the
quadratic twist square class without enumerating every nonzero scalar multiple.
With `--monic-only`, only leading coefficient `1` is used. Before point
counting, the search applies presentation reduction. This skips many repeated
presentations, but it is not a full isomorphism-class computation.

Before generating a full reduction orbit, the search checks whether the
presentation is already equivalent to an accepted curve, then applies the
Hasse-Witt prefilter. This avoids PGL2 orbit generation for many candidates that
already cannot have a trinomial L-polynomial. Point counts are computed over
`F_{p^k}` for `k = 1, ..., g` only after those filters pass, and the
L-polynomial coefficients are recovered from Newton identities. The
implementation rejects a candidate immediately when one of the required zero
coefficients `a_1, ..., a_{g-1}` is nonzero, so it does not count points over
larger extensions for candidates that already cannot have a trinomial
L-polynomial.

Both implementations cache point-counting data for each extension field when it
is first needed: the field model, square-count table, and powers `x^e` for all
field elements and exponents up to `2g + 2`. Candidate evaluation then sums
only the nonzero terms of `f(x)` using those cached power tables.

The C implementation also uses FLINT for selected polynomial operations where
the general-purpose library is expected to help: squarefree checks via
derivative/gcd over `F_p`, irreducible modulus construction for extension
fields, and nontrivial Hasse-Witt polynomial powers.
