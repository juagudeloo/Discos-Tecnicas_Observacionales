# PhotUN 

---

This code was developed to apply photometric measurements on the channels 1,2,3,4 for .fits files of the IRAC instrument and on the W1,W2,W3,W4 channels for .fits files of the WISE instrument.

## photUN_SPITZER.py

This script deals with SPITZER channels

### SPITZER_notebook.ipynb

This jupyter notebook develops and explains all the steps done to obtain the photUN_SPITZER.py code to work.

## photUN_WISE.py

Still in development...

### SPITZER_notebook.ipynb

Still in development...

## General parameters

The scripts receive the following parameters (you can check them also by typing `python3 photUN_SPITZER.py --help` or `python3 photUN_WISE.py --help`).

```
positional arguments:
  {IRAC,WISE}           String type that specifies 'IRAC' or 'WISE' for the instrument to analyze.

options:
  -h, --help            show this help message and exit
  -fd FITS_DIR, --fits_dir FITS_DIR
                        A string type that indicates the path to the .fits files.
  -cad CATALOG_DIR, --catalog_dir CATALOG_DIR
                        A string type that indicates the path to the GAIA DR2 catalog.
  -od OUTPUT_DIR, --output_dir OUTPUT_DIR
                        A string type that indicates the path to save the resulting .csv output file.
  -r R, --circular_aperture R
                        Integer value that determines the circular aperture in pixels.
  -r_in R_IN, --annular_int_radius R_IN
                        Integer value that determines the annular aperture internal radius in pixels.
  -r_out R_OUT, --annular_ext_radius R_OUT
                        Integer value that determines the annular aperture external radius in pixels.
  -v, --verbose

It is intended to work with MIPS in the future.
```