# Prediction of viral fitness evolution

## Summary
Since our model simply represents the viral lineage growth parameter ($\beta_l$) as the linear combination of the effects of S substitutions, the model can predict the total effect of a set of substations on relative Re. Using this property of the model, we predicted the ancestral relative viral fitness for each node in the BA.5 tree. We calculated these values as the sum of the posterior means of the effects of substitutions of interest. To reconstruct the ancestral relative viral fitness of each node of the BA.5 tree, we first reconstructed the ancestral state of the S substitution profile in each node of the tree using a parsimony method, implemented by the castor package. Subsequently, we predicted the relative viral fitness for each node according to the reconstructed ancestral mutation profile for the node.



