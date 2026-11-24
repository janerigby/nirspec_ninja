# NIRSpec finder charts

`nirspec_finder` overlays the four MSA quadrants, IFU and fixed slits on an
image, using the pointing and orientation of a NIRSpec exposure. It helps
identify sources and sky within each aperture. The display aligns the mean
MSA edges with the plot axes and places the IFU to the left.

## Inputs

The positional input must be a **NIRSpec exposure FITS file** with pointing
metadata, such as a RATE file. A NIRCam image belongs in `--image`, not in the
positional input. Required header values are either
`V2_REF/V3_REF/RA_REF/DEC_REF/ROLL_REF` or `RA_V1/DEC_V1/PA_V3`; SCI values take
precedence over primary-header values.

With `--image`, supply a 2D celestial-WCS FITS image, such as a NIRCam mosaic.
Partial image coverage is supported: uncovered areas are blank, and aperture
outlines remain visible wherever they intersect the chart. At least some
finite image data must fall within the displayed field.

Without `--image`, the tool retrieves WISE W1 (3.4 μm) imagery from IRSA.
This needs internet access on the first retrieval; downloaded cutouts are cached.

## Examples

Full NIRSpec field over a local image:

```sh
python -m nirspec_ninja.finder.nirspec_finder \
  exposure_rate.fits \
  --image nircam_mosaic.fits \
  --output finder.png --regions
```

Zoom around the fixed slits using the same files:

```sh
python -m nirspec_ninja.finder.nirspec_finder \
  exposure_rate.fits \
  --image nircam_mosaic.fits \
  --output finder_FSzoom.png --regions --FSzoom
```

Use WISE rather than a local image:

```sh
python -m nirspec_ninja.finder.nirspec_finder exposure_rate.fits \
  --output wise_finder.png --regions
```

## Options and products

| Option | Purpose |
| --- | --- |
| `--image PATH` | Use a local celestial-WCS FITS image instead of WISE. |
| `--output PATH` | Save the PNG here; default is `<input_stem>_finder.png` beside the exposure. |
| `--regions` | Also write DS9 aperture polygons to a `.reg` file. |
| `--FSzoom` | Show approximately 100 × 60 arcseconds around the IFU and S200A1, S200A2, S400A1 and S1600A1. |
| `--size VALUE` | Set field width in arcminutes. Full view is square; FS zoom retains a 5:3 width-to-height ratio. |
| `--cache PATH` | Choose the WISE cache directory; default is `finder_cache` beside the package's finder module. |

MSA outlines are cyan, IFU magenta, and fixed slits yellow. The chart includes
an N/E compass, field size and S200A1 center coordinates. A companion `.json`
file records the input, image source, pointing, SIAF version, orientation and
field dimensions. Output files sharing the same basename are overwritten;
use different `--output` names for full and zoom views.

S200B1 lies elsewhere in the field and is outside the default FS zoom.
Use the full view to see it. Region files include all ten apertures even when
some lie outside the displayed field.

## Warnings and limitations

- A `pysiaf` PRD version warning means the installed aperture-reference tables
  differ from the latest release. It does not itself stop chart generation.
- Astropy `datfix` and `obsfix` warnings describe conversions of date and
  observatory-position header values; these do not themselves cause a failure.
- The tool uses one exposure's pointing. It does not reconstruct constituent
  pointings or dithers from a combined product.
- Full-field charts from fine-scale NIRCam images can require substantial
  memory. `--FSzoom` reduces the displayed area and the corresponding image grid.
