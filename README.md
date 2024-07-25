# PhotUN 

---

This code was developed to apply photometric measurements on the channels 1,2,3,4 for .fits files of the IRAC instrument and on the W1,W2,W3,W4 channels for .fits files of the WISE instrument.

## photUN.py

This script deals with IRAC or WISE channels.

### SPITZER_notebook.ipynb

This jupyter notebook develops and explains all the steps done to obtain the `photUN.py IRAC` code to work.

### WISE_notebook.ipynb

This jupyter notebook develops and explains all the steps done to obtain the `photUN.py WISE` code to work.

## General parameters

The script receive the following parameters (you can check it also by typing `python3 photUN.py --help`).

```
usage: photUN.py [-h] -fd FITS_DIR -cad CATALOG_DIR [-od OUTPUT_DIR] [-r R]
                 [-r_in R_IN] [-r_out R_OUT] [-v]
                 {IRAC,WISE}

process IRAC (SPITZER) or WISE .fits files to get their photometry.

positional arguments:
  {IRAC,WISE}           String type that specifies 'IRAC' or 'WISE' for the
                        instrument to analyze.

options:
  -h, --help            show this help message and exit
  -fd FITS_DIR, --fits_dir FITS_DIR
                        A string type that indicates the path to the .fits
                        files.
  -cad CATALOG_DIR, --catalog_dir CATALOG_DIR
                        A string type that indicates the path to the GAIA DR2
                        catalog.
  -od OUTPUT_DIR, --output_dir OUTPUT_DIR
                        A string type that indicates the path to save the
                        resulting .csv output file.
  -r R, --circular_aperture R
                        Integer value that determines the circular aperture in
                        pixels.
  -r_in R_IN, --annular_int_radius R_IN
                        Integer value that determines the annular aperture
                        internal radius in pixels.
  -r_out R_OUT, --annular_ext_radius R_OUT
                        Integer value that determines the annular aperture
                        external radius in pixels.
  -v, --verbose

It is intended to work with MIPS in the future.
```