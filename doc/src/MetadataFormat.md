# MetadataFormat

## Overview

`metadata.csv` is used by `buildIS` in training mode, while `metadata_test.csv` is used by `exploreIS` in evaluation mode. Both files share an identical column layout.

## Column Schema

### Instances
The first column contains the instance label/identifier (string), with one row per instance.

### Feature Columns (`feature_*`)
Numeric columns representing instance features. In the reference dataset, there are 10 such columns.

### Algorithm Columns (`algo_*`)
Numeric columns representing an algorithm's raw performance value on that instance. In the reference dataset, there are 10 algorithm columns (NB, LDA, QDA, CART, J48, KNN, L_SVM, poly_SVM, RBF_SVM, RandF).

### Optional Source Column
An optional `source` column is used for per-instance provenance labels (not present in the reference dataset).

## Configuration and Constraints

### Column Selection
`opts.selvars.feats` and `opts.selvars.algos` in `options.json` can restrict used columns by name.

### Evaluation Mode Requirements
For `metadata_test.csv` in evaluation mode:
- Feature columns must match the training `metadata.csv` in both set and order (validated by `InstanceSpace.explore()/evaluateTestSet`).
- Algorithm columns: Known algorithms match their trained column position; new algorithms present in the test file are appended as new columns.

## Example

Confirmed via direct inspection of test/data/metadata.csv and test/data/metadata_test.csv's header row (identical in both files):

```
Instances,feature_Max_Normalized_Entropy_attributes,feature_Normalized_Entropy_Class_Attribute,feature_Mean_Mutual_Information_Attribute_Class,feature_ErrorRate_Decision_Node,feature_WeightedDist_StdDev,feature_Max_Feature_Efficiency_F3,feature_Collective_Feature_Efficiency_F4,feature_Training_Error_Linear_Classifier_L2,feature_Fraction_Points_Class_Boundary_N1,feature_Nonlinearity_Nearest_Neighbor_Classifier_N4,algo_NB,algo_LDA,algo_QDA,algo_CART,algo_J48,algo_KNN,algo_L_SVM,algo_poly_SVM,algo_RBF_SVM,algo_RandF
```

First data row (metadata.csv), showing the value format:

```
abalone,0.332548382,0.131686939,0.229928422,0.386458668,2.343235444,0.076929217,0.095203899,0.260234614,0.404836009,0.327388088,0.280095763,0.239715007,0.25401222,0.24533949,0.256491398,0.256969153,0.230584159,0.277410974,0.221197935,0.234432252
```
