# Hyperelliptic Curve Finder

Find presentations of hyperelliptic curves

```text
y^2 = f(x)
```

over odd prime fields `F_p` whose L-polynomial is a trinomial

```text
L(t) = 1 + a_g t^g + p^g t^{2g}
```

with nonzero middle coefficient `a_g`.

The program lists reduced presentations, not certified isomorphism classes.

## Usage

```bash
python3 hyperelliptic_finder.py p g [--max N] [--output FILE] [--reduction {pgl2,affine}] [--quiet]
```

Arguments:

- `p`: odd prime characteristic. Characteristic `2` is not supported for models `y^2 = f(x)`.
- `g`: genus.
- `--max N`: maximum number of presentations to save. Use `0` for the complete search.
- `--output FILE`: output path. Defaults to `hyperelliptic_results.txt`.
- `--reduction {pgl2,affine}`: presentation reduction before point counting. Defaults to `pgl2`.
- `--quiet`: suppress per-presentation progress messages.

If `--output` ends in `.txt`, a matching `.json` file is also written. If it
ends in `.json`, a matching `.txt` file is also written.
Progress printed to the terminal is also saved to a matching `.log.txt` file.

## Examples

```bash
python3 hyperelliptic_finder.py 3 1 --max 10
python3 hyperelliptic_finder.py 3 2 --output genus2.txt
python3 hyperelliptic_finder.py 5 2 --max 0 --output complete-results.txt
python3 hyperelliptic_finder.py 3 2 --reduction affine --output affine-results.txt
```

## Batch Runs

Use `test.py` to run many searches and save each case to a separate file.
By default it runs primes `3, 5, 7` and genera `1` through `5`.

```bash
python3 test.py --outdir batch_results --resume --timeout 600
```

Useful options:

- `--outdir DIR`: directory for per-case `.txt`/`.json` output files.
- `--reduction {pgl2,affine}`: reduction mode passed to `hyperelliptic_finder.py`.
- `--max N`: maximum curves per case; `0` means complete search.
- `--timeout SECONDS`: maximum runtime per case; `0` means no timeout.
- `--min-genus G` and `--max-genus G`: restrict the genus range.
- `--resume`: skip cases whose JSON output already exists.
- `--quiet`: pass `--quiet` to `hyperelliptic_finder.py`.

The script also writes:

```text
batch_results/batch_summary.json
batch_results/batch.log.txt
```

`batch_summary.json` is updated after each completed case, so partial progress
is preserved if the batch is interrupted with Ctrl+C. Each individual
`hyperelliptic_finder.py` run also saves its own results and progress log
incrementally.

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

For batch results:

```bash
python3 lmfdb.py batch_results/*.json --out lmfdb_accuracy.json --delay 1.0 --timeout 30
```

Useful options:

- `--out FILE`: output JSON report.
- `--delay SECONDS`: wait between LMFDB requests.
- `--timeout SECONDS`: HTTP timeout per request.

The comparison is by L-polynomial/isogeny class, not by individual curve
presentation. LMFDB may not have records for every genus and prime searched.
The reported `accuracy_among_available` only uses LMFDB records that were found.

## Output

Detailed results are saved to the text and JSON output files. The terminal only
prints run status, per-presentation progress, reduction skips, early rejections,
and saved matches. Use `--quiet` to print only saved matches and the final
summary. The same terminal progress is written to the `.log.txt` file.

Each saved curve presentation is indexed and includes:

- `f(x)` for the presentation `y^2 = f(x)`.
- The coefficient list for `f(x)`, in ascending order `[c_0, c_1, ...]`.
- The middle L-polynomial coefficient `a_g`.

Only the middle coefficient of the L-polynomial is shown because the other
nonzero coefficients in the trinomial form are fixed as `1` and `p^g`.

Progress is saved after every matching presentation. If the program is
interrupted with Ctrl+C, all matches found before the interruption remain in the
output files.

## Reduction Options

The `--reduction` option controls which repeated presentations are skipped
before point counting.

### `--reduction pgl2`

This is the default. It uses fractional linear transformations

```text
x -> (ax + b)/(cx + d)
```

with `ad - bc != 0` over `F_p`. This includes affine transformations and also
transformations such as inversion `x -> 1/x`. It can identify more equivalent
presentations, including some cases where a degree `2g + 1` presentation and a
degree `2g + 2` presentation are related.

### `--reduction affine`

This uses only transformations

```text
x -> ax + b
```

with `a != 0` over `F_p`. It is cheaper than `pgl2`, but removes fewer repeated
presentations.

Both modes normalize transformed equations back to monic form only when the
corresponding `y`-scaling exists over `F_p`.

## Notes

The search enumerates monic squarefree polynomials of degree `2g + 1` and
`2g + 2` over `F_p`. Before point counting, it applies presentation reduction.
This skips many repeated presentations, but it is not a full isomorphism-class
computation.

Point counts are computed over `F_{p^k}` for `k = 1, ..., g`, and the
L-polynomial coefficients are recovered from Newton identities before checking
the trinomial condition. The implementation rejects a candidate immediately when
one of the required zero coefficients `a_1, ..., a_{g-1}` is nonzero, so it does
not count points over larger extensions for candidates that already cannot have
a trinomial L-polynomial.
