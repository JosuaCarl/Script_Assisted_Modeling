# Script_Assisted_Modeling
Python scripts, that aid in the construction of genome-scale metabolic models (GEMs), after a draft with CarveMe.

# Python scripts
- dependencies on COBRAPy, libSBML, pandas, requests, os, memote, numpy, tqdm, bioservices, gffpandas
  - most are available via pip

# Databases
Some scripts require the BiGG, MetaNetX and SEED Databases, structured like this, starting form script location:
- ../Databases/BiGG
  - bigg_models_metabolites.tsv
  - bigg_models_reactions.tsv
  - accessible on http://bigg.ucsd.edu/data_access
- ../Databases/MetaNetX
  - chem_prop.tsv
  - reac_prop.tsv
  - accessible on https://www.metanetx.org/mnxdoc/mnxref.html
- ../Databases/SEED
  - compounds.tsv
  - reactions.tsv
  - accessible on https://github.com/ModelSEED/ModelSEEDDatabase/tree/master/Biochemistry

  ```
  mkdir ../Databases ../Databases/BiGG ../Databases/MetaNetX ../Databases/SEED
  cd ../Databases/BiGG
  wget http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt
  wget http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt
  mv bigg_models_reactions.txt bigg_models_reactions.tsv
  mv bigg_models_metabolites.txt bigg_models_metabolites.tsv
  cd ../MetaNetX
  wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv
  wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_prop.tsv
  cd ../SEED
  wget https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/reactions.tsv
  wget https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/compounds.tsv
  ```

# iPython notebooks
- made with jupyter-lab
