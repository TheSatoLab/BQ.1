# Prediction of viral fitness evolution

## Summary
Since our model simply represents the viral lineage growth parameter ($\beta_l$) as the linear combination of the effects of S substitutions, the model can predict the total effect of a set of substations on relative Re. Using this property of the model, we predicted the ancestral relative viral fitness for each node in the BA.5 tree. We calculated these values as the sum of the posterior means of the effects of substitutions of interest. To reconstruct the ancestral relative viral fitness of each node of the BA.5 tree, we first reconstructed the ancestral state of the S substitution profile in each node of the tree using a parsimony method, implemented by the castor package. Subsequently, we predicted the relative viral fitness for each node according to the reconstructed ancestral mutation profile for the node.

## Note
Please download the GISAID Metadata from the sction "Download packages" (https://gisaid.org/) and save it in this directory.


## Usage
```bash

cd <THIS DIRECTORY>

#extract mutation info

python3 script/summarize_mut_info.ver2.py \
  metadata.txt \
  metadata.mut_long.txt

#predict viral fitness evolution

R --vanilla --slave --args \
  metadata.txt \
  metadata.mut_long.txt \
  < script/predict_viral_fitness_evolution.R

```

## System requirements
* **R** v4.1.2
* **tidyverse** v1.3.1
* **data.table** v1.14.2
* **ggplot2** v3.3.6
* **ape** v5.6.2
* **ggtree**  v3.2.1
* **castor** v1.7.2
* **phangorn** v2.9.0
* **patchwork** v1.1.1
* **RColorBrewer** v1.1.3
* **ggnewscale** v0.4.7
* **scales** v1.2.0
