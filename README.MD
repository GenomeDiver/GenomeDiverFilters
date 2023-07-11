Phenotypic filters for Genome Diver. 
The executable is ran by Genome Diver via Nextflow. 

### Usage 
- - - 
Run sbt packageBin to create a jar file from sources.  
Run run-blank-example.sh to see parameters 

### Dependencies
- - - 
Filters depend on Exomiser output (TSV) as well the Lirical data folder. 

### Discussion
- - - 
Genes/Diseases isolated by Exomiser are associated with phenotypes. These phenotypes first filtered by the input. Afterwards they are ranked following via Lirical in terms of likelihood ratios to their respective diseases. Phenotypes are further ranked by their uniqueness to a particular disease via dispersion index. 

The resulting phenotypes satisfy the criteria of being unique to the input as well as discerning for downstream processing should they be accepted by the caregiver. 

## Thanks
Nathan Pearson (initial proposals on filtering criteria)  
Faygel Beren (refinement and incorporation of HPO frequency) 
