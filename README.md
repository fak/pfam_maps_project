This module combines a script to project a list of validated domains described in the [InCoB 2012 conference supplement](http://www.biomedcentral.com/bmcbioinformatics/supplements) onto the target dictionary of a given ChEMBL release. It also includes a facility to manually curate conflicting mappings through a Django curation interface and feed these corrections back into the database. 

This is to replace an older mechanism that generated mappings for each release of the ChEMBL database but needed to be adjusted with the schema change that occurred at the transition between chembl_14 and  chembl_15. Mappings of older versions are provided on a dedicated [website](http://www.ebi.ac.uk/~fkrueger/mapChEMBLPfam/).

