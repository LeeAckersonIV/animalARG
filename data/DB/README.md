# This directory contains the reference databases used for ARG and 16s rRNA alignement
There are two default databases utilized by this pipeline, both of which have already been indexed via BWA.
Additionally, both databases have been clustered individually via cd-hit-est to remove redundancy.

## Utilized reference databases:
### The Comprehensive Anitbiotic Resistance Database (CARD)
To utilize your own reference and associated index files:
 ```sh
  nextflow run main.nf --CARD_db 'your/path/to/db_16SrRNA.files'
  ```
### GreenGenes 16s rRNA Database
To utilize your own reference and associated index files:
 ```sh
  nextflow run main.nf --GG_db 'your/path/to/db_ARG.files'
  ```
