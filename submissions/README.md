# Submssions
This directory contains the scripts and data to prepare assembled consensus genomes for submission to GISAID, GenBank, and WA DOH.

## Data

Files stored in the `source-data/` directory are:

* `authors/` subdirectory contains separate files that lists all the authors for each sample source.
These files are named with their respective sources:
    * `altius.txt`
    * `scan.txt`
    * `sfs.txt`
    * `wa-doh.txt`

* `us_state_abbreviations.csv` contains US states and their abbreviations.
This is used to create the appropriate strain name based on the sample's original state.

* `variants_of_concern.tsv` contains mapping of variants of concerns for the clade name from Nextclade, WHO, and Pangolin.

* `washington_counties.txt` contains a list of all counties in the state of Washington.
This will be used to flag samples with counties that are from outside of Washington.

* `washington_pumas_to_county.csv` contains a mapping of Washington PUMAs to their counties.
This file only contains PUMAs that map to **one** specific county within Washington.
This is taken from the PUMA maps listed at
https://www.census.gov/geographies/reference-maps/2010/geo/2010-pumas/washington.html



