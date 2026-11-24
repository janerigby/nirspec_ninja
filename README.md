# nirspec_ninja

## Overview

**`nirspec_ninja` is a set of tools for JWST NIRSpec data reduction**. Things the ninja can do:

- [Query and download large data sets](src/nirspec_ninja/archive/README.md), with retries to handle timeouts.
- [Process large amounts of NIRSpec data on multiple CPUs](src/nirspec_ninja/pipeline/README.md).
- [Reduce serendipitous NIRSpec fixed-slit spectra](src/nirspec_ninja/fixed_slit/README.md). The fixed slits are always open; while standard MSA and IFU reduction ignores them, the fixed slits contain useful information about sky brightness.
- [Estimate and subtract sky backgrounds](src/nirspec_ninja/subtract_background/README.md), including using the fixed slits to estimate sky
- [Detect artifacts](src/nirspec_ninja/artifacts/README.md), specifically glow from MSA electrical shorts and column artifacts from problems with 1/f noise correction.
- [Make custom masks of illuminated pixels](src/nirspec_ninja/masks/README.md) for 1/f noise correction. Mask values of True/1 identify usable background pixels; False/0 identifies excluded pixels. See [validation coverage](src/nirspec_ninja/masks/VALIDATION_COVERAGE.md) for the grating/filter combinations and detectors checked so far.
- [Measure quality of data reduction](src/nirspec_ninja/diagnostics/README.md), including impact of 1/f correction and presence of gradients in fixed-slit spectra.
- [Make finder charts showing NIRSpec apertures](src/nirspec_ninja/finder/README.md), using WISE or local images.
- [Patch bugs in the JWST pipeline's 1/f noise correction (present in v3.0.0)](src/nirspec_ninja/pipeline/patches/README.md). Importing the package does not apply these patches; follow the linked instructions to apply them in your reduction environment.

## Contributors

Individuals: Jane Rigby, Brian Welch, Rosalia O'Brien, and Taylor Hutchison  
Teams: TEMPLATES, LEGGOS

## Getting started

See [INSTALL.md](INSTALL.md) for installation, package imports, and notebook setup.
The links above provide usage instructions for each tool.

## Development notes

See [migration status and development notes](MIGRATION_STATUS.md) for implementation history, compatibility changes, and remaining work, and [live-service test instructions](scripts/README.md) for WISE and AWS/MAST smoke tests.

## GAI acknowledgement: 
Much of this code was written or adapted using the Codex code assistant version 6.1
