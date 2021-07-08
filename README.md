# PTM-Site-Validator
PTM site validator: Amino acid sites found to have post-translational modification (PTM) from Proteome Discoverer searches can be validated by identifying modification site specific immonium ion. Note: The tool requires raw files used for the search in .mgf format in the same folder where the PTM site validator is present. 

For any enquiries please contact:
sandeep.kolya@gmail.com or chinnu.kemmaai@gmail.com

How to use PTM-Site-Validator

Command line argument:
```
ptm_site_validator> python PTM_site_validator.py -i PSMs.txt -ft 0.02 -r Path_to_raw_data

```
positional arguments:
  -i          Path to the PSM file from PTM search in Proteome Discoverer.
  -ft         Maximum search m/z range for the identification of immonium ion
              in MS/MS spectra. It is recommended to keep the same search
              tolerance used in the main search.
  -r          Path to the folder in which raw files are stored in .mgf format.

optional arguments:
  -h, --help  show this help message and exit
```
