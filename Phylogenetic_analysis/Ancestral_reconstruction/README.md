# Inferring branches where substitution events occurred

## Summary
To infer the branches where substitution events occurred at the five convergent sites (346, 444, 452, 460, 486) in the trees of Omicron lineages, we reconstructed the ancestral state of the substitution profile at the convergent sites in each node using a parsimony method, implemented by the castor package (https://www.rdocumentation.org/packages/castor/versions/1.7.2). Internal nodes with substitution probabilities above or equal to 0.5 were annotated as having the substitution. Branches where substitutions took place for each site were denoted as branches connecting an ancestral internal node with no substitution to an internal node that has a substitution. Additionally, 70% of tips descending from that internal node were also required to have the substitution and at least 3 tips needed to be descended from the node, to avoid picking up branches with low support or clades that reverted back to the original residue.

## Note
GISAID Metadata downloaded from the sction "Download packages" is needed as an input.

## Usage
```bash

cd <THIS DIRECTORY>

#extract mutation info

python3 script/summarize_mut_info.ver2.py \
  metadata.txt \
  metadata.mut_long.txt

#Inferring branches where substitution events occurred and visualizing it
R --vanilla --slave --args \
  metadata.txt \
  metadata.mut_long.txt \
  < script/detect_transition_branch_and_visualization.R

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




