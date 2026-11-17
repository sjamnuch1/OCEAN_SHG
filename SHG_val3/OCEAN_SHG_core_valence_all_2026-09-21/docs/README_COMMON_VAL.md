# OCEAN Common-only valence input processor

`ocean_common_val.pl` is a restricted adaptation of OCEAN's `ocean.pl`.
It performs only the Common stage and requires the processed calculation
mode to be `val`. By default it creates `Common_VAL`, preserving an existing
core `Common` directory.

## Usage

From the run root:

```bash
perl /path/to/ocean_common_val.pl \
  --input valence.in \
  --root "$PWD" \
  --ocean-bin /path/to/OCEAN/bin
```

`OCEAN_BIN` can be used instead of `--ocean-bin`. The output defaults to
`ROOT/Common_VAL`. `--output` may select another fresh direct child of the
run root. Existing files or directories are never overwritten.

The input must contain the normal OCEAN setting:

```text
calc.mode { val }
```

For the finite-q SHG valence source, it must also specify the required small
q using the syntax supported by your OCEAN version. The wrapper warns when
the resulting `calc.photon_q` is zero because normal q=0 valence PREP is
valid, but `shg_xc_minus_valence.x` divides the overlap by q.

The wrapper runs the same Common helpers as the inspected original launcher:

1. `parseJSON.pl` with `oparse.json` and `oparse.type.json`, falling back to
   legacy `parse.pl`/`oparse.h` when parseJSON is unavailable.
2. `extractPsp.pl` and `psp_parser.pl`.
3. `defaults.pl`, `structure.pl`, `edges.pl`, and `photon.pl`.
4. Validation of `postDefaultsOceanDatafile`, including `calc.mode=val`.

It then writes `common_val.complete`. It does not create or launch OPF, DFT,
PREP, SCREEN, CNBSE, or zWFN. As in the original launcher, failure of
`psp_parser.pl` is a warning; all other helper failures are fatal.

## Connection to the valence workflow

The next stages are explicit:

```bash
# Generate the matching finite-q DFT data using the valence Common settings.
# Then perform the separate PREP conversion:
perl prep_val.pl \
  --root /path/to/core_run \
  --common /path/to/core_run/Common_VAL/postDefaultsOceanDatafile \
  --dft-dir /path/to/valence_DFT \
  --ocean-bin /path/to/OCEAN/bin
```

The Common-only wrapper does not generate shifted DFT wavefunctions. The DFT
run supplied later must match `kmesh`, `kshift`, `photon_q`, bands, structure,
pseudopotentials, and spin settings from `Common_VAL`.

The earlier `prepare_valence.pl` also needs its source paths changed from
`Common`, `DFT`, and `PREP` to the selected valence locations before consuming
this side-by-side workflow.

## Verification

```bash
perl -c ocean_common_val.pl
python3 test_ocean_common_val.py
```

The mocked integration test verifies Common-only execution, val-mode
enforcement, overwrite protection, the completion marker, and absence of
later-stage directories. It does not replace a test with the real OCEAN
parser and a production input file.

The upstream OCEAN license is included. This wrapper invokes installed OCEAN
helpers; it does not reimplement their input/default logic.
