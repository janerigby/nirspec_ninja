# Installation and use during migration

Activate the Python environment used for your reductions, then install once
in each environment:

```sh
python -m pip install --no-deps --no-build-isolation -e /Users/jrrigby1/AI_experiments/Coding/nirspec_ninja
```

The `-e` installation uses the code directly from this directory, so code edits
are available without reinstalling. It is already installed in
`/Users/jrrigby1/AI_experiments/References_and_Tools/JWSTDP-3.0.0-10-py314-macos-arm64`.
The command uses your existing astronomy dependencies without upgrading them.

Import the migrated functions:

```python
from nirspec_ninja.masks import make_custom_NIRSpec_mask
from nirspec_ninja.fixed_slit.source_detection import wrap_detect_source_s2d
from nirspec_ninja.pipeline.patches.nirspec_flicker_patch import install_patch
from nirspec_ninja.pipeline import notebook_workers
from nirspec_ninja.pipeline.associations import writel3asn
```

Use the package imports above. Known Coding callers have been updated and
unused import-only redirects retired. Old command-line entry points remain
available where retained.

For notebooks, select a kernel using the same Python environment. Restart the
kernel, Python session, or running workers after code changes so they load the
updated modules. Importing the patch module does NOT apply the patch;
calling `install_patch()` applies it to that running Python process.

To run the package regression tests, from the nirspec_ninja directory:

```sh
python -m unittest discover -s tests -v
```
