# prep_val.pl

Adaptation of the supplied original OCEAN PREP/prep.pl (2026-08-19 source
archive). It runs the normal valence PREP conversion into a fresh
`ROOT/PREP_VAL/BSE`, with its own `ROOT/PREP_VAL/prep.json`.
The existing `ROOT/PREP` and upstream files are not modified.

## Run

If Common and DFT already describe the completed valence DFT calculation:

```bash
perl /path/to/prep_val.pl --root /path/to/run --ocean-bin /path/to/OCEAN/bin
```

To keep the original core Common/DFT directories unchanged, supply an
independently prepared valence configuration and matching DFT directory:

```bash
perl /path/to/prep_val.pl \
  --root /path/to/core_run \
  --common /path/to/valence_run/Common/postDefaultsOceanDatafile \
  --dft-dir /path/to/valence_run/DFT \
  --core-common /path/to/core_run/Common/postDefaultsOceanDatafile \
  --edge-index 1 \
  --ocean-bin /path/to/OCEAN/bin \
  --dry-run
```

Remove `--dry-run` to launch `ocean_prep.x`. Alternatively, the selected
inputs can live in `core_run/Common_VAL` and `core_run/DFT_VAL`.
The installed OCEAN bin directory must contain ocean_prep.x, specpnt.5,
Pquadrature, and sphpts. OCEAN_BIN is used if --ocean-bin is omitted;
otherwise the directory containing prep_val.pl is the fallback.
The script is invoked from the run root (or with --root), not inside PREP.

The selected Common JSON must have `calc.mode=val`. The script deliberately
does not change a core configuration's mode or pretend that its DFT has
the missing occupied shifted states. DFT's BSE stage must be complete;
kmesh, kshift, photon_q and nonzero_q must match between Common and DFT.
The original DFT directory naming and split/combined QE logic are retained.
ABINIT non-split input follows the original path; split ABINIT is rejected
because the inherited link routine cannot handle it. Exact JSON metadata
matching is conservative and may reject numerically equivalent encodings.

## What it does

- Requires val mode but enables both `prep.tmels` and `prep.cks`. This produces
  the finite-q valence source data and same-gauge core projectors in one PREP.
- Imports one edge from `--core-common` using the 1-based `--edge-index`
  (default 1), or accepts an explicit `--edge 'SITE N L'`.
- Writes the usual PREP inputs and links the matching DFT wavefunctions read-only.
- Runs ocean_prep.x using the configured computer.para_prefix launcher.
- Checks the expected output files exist and are nonempty before recording completion.
- Stores selected Common and DFT JSON snapshots in PREP_VAL.
- Refuses an existing PREP_VAL, including symlinks; there is no overwrite option.

The expected BSE files include val.u2.dat, con.u2.dat, ptmels.dat,
tmels.info, enkfile, wvfcninfo, wvfvainfo, parcksv.*, and parcksc.*. A failure leaves its partial
PREP_VAL and logs for inspection; move that directory aside before retrying.
Support for restart/cache reuse was intentionally disabled to keep a
previous run or a core directory from being overwritten.

## What it does not do

PREP does not perform DFT. For a finite-q source, generate the matching
shifted occupied/conduction wavefunctions upstream first. A q=0 run is
allowed as in normal PREP, but the earlier shg_xc_minus_valence driver
requires a nonzero small q for its dipole construction.

`prepare_valence.pl` consumes this output with `--prep-dir PREP_VAL`; no
renaming or replacement of the core PREP tree is needed. It stages the
projectors and generates all three Cartesian `mels` files immediately.

## Verification

```bash
perl -c prep_val.pl
python3 test_prep_val.py
```

Perl syntax and a mocked end-to-end test passed. The test checks the val-only
guard, q mismatch, dry-run behavior, split QE links, generated VAL flags,
completion checks, and preservation of core/upstream data. It does not
validate numerical wavefunctions: real ocean_prep.x execution requires
your installed OCEAN build and DFT dataset.

The upstream OCEAN license is included. Existing PREP conversion and file
layout logic are retained rather than replacing the numerical algorithms.
