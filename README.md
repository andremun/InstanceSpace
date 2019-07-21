# Instance Space Analysis: A toolkit for the assessment of algorithmic power

Instance Space Analysis is a methodology for the assessment of the strengths and weaknesses of an algorithm, and an approach to objectively compare algorithmic power without bias introduced by restricted choice of test instances. At its core is the modelling of the relationship between structural properties of an instance and the performance of a group of algorithms. Instance Space Analysis allows the construction of **footprints** for each algorithm, defined as regions in the instance space where we statistically infer good performance. Other insights that can be gathered from Instance Space Analysis include:

-	Objective metrics of each algorithm’s footprint across the instance space as a measure of algorithmic power;
-	Explanation through visualisation of how instance features correlate with algorithm performance in various regions of the instance space;
-	Visualisation of the distribution and diversity of existing benchmark and real-world instances;
-	Assessment of the adequacy of the features used to characterise an instance;
-	Partitioning of the instance space into recommended regions for automated algorithm selection;
-	Distinguishing areas of the instance space where it may be useful to generate additional instances to gain further insights.

The unique advantage of visualizing algorithm performance in the instance space, rather than as a small set of summary statistics averaged across a selected collection of instances, is the nuanced analysis that becomes possible to explain strengths and weaknesses and examine interesting variations in performance that may be hidden by tables of summary statistics.

This repository provides a set of MATLAB tools to carry out a complete Instance Space Analysis in an automated pipeline. It is also the computational engine that powers the Melbourne Algorithm Test Instance Library with Data Analytics ([MATILDA](http://matilda.unimelb.edu.au/)) web tools for online analysis. For further information on the Instance Space Analysis methodology can be found [here](http://matilda.unimelb.edu.au/our-methodology).

## Installation Instructions

The only requirement for the software to run is to have a current version of [MATLAB](http://www.mathworks.com). It has been tested and known to work properly in version r2018b. Earlier versions may fail to support the processing of ```.json``` files. To make use of all functionalities, the [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) library must be installed. This repository contains the Windows binaries of LIBSVM, but other platforms require compilation from source code. We suggest adding the software folder to the MATLAB path.

## Working with the code

The main interfase is the script ```example.m``` which provides the path for the ```metadata.csv``` file, and constructs the ```options.json``` file. The path provided will also be the location of all the software outputs, such as images (in ```.png``` format), tables (in ```.csv``` format) and raw intermediate data (in ```.mat``` format).

## The metadata file

The ```metadata.csv``` file should contain a table where each row corresponds to a problem instance, and each column must strictly follow the naming convention mentioned below:

-	**instances** instance identifier - We expect instance identifier to be of type "String". This column is mandatory.
-	**source** instance source - This column is optional
-	**feature_name** The keyword "feature_" concatenated with feature name. For instance, if feature name is "density", header name should be mentioned as "feature_density". If name consists of more than one word, each word should be separated by "_" (spaces are not allowed). There must be more than two features for the software to work. We expect the features to be of the type "Double".
-	**algo_name** The keyword "algo_" concatenated with algorithm name. For instance, if algorithm name is "Greedy", column header should be "algo_greedy". If name consists of more than one word, each word should be separated by "_" (spaces are not allowed). You can add the performance of more than one algorithm in the same ```.csv```. We expect the algorithm performance to be of the type "Double".

Moreover, empty cells, NaN or null values are not allowed. We expect you to handle missing values in your data before processing. You may use [this file](https://matilda.unimelb.edu.au/matildadata/graph_coloring_problem/metadata.csv) as reference.

## Options

The script ```example.m``` contains code that prepares all the settings used by the code. Broadly, there are settings required for the analysis itself, settings for the pre-processing of the data, and output settings. For the first these are divided into general, dimensionality reduction, algorithm selection and footprint construction settings. For the second, the toolkit has routines for bounding outliers, scale the data and select features.

### General settings

-	```opts.perf.MaxMin``` determines whether the algorithm performance values provided are **efficiency** measures that should be maximised (set as ```TRUE```), or **cost** measures that should be minimised (set as ```FALSE```).
-	```opts.perf.AbsPerf``` determines whether good performance is defined absolutely, e.g., misclassification error is lower than a 20%, (set as ```TRUE```), or if it is defined relatively to the best performing algorithm, e.g., misclassification error is within at least 5% of the best algorithm, (set as ```FALSE```).
-	```opts.perf.epsilon``` corresponds to the threshold used to calculate good performance. It must be of the type "Double".
-	```opts.general.betaThreshold``` corresponds to the fraction of algorithms in the portfolio that must have good performance in the instance, for it to be considered an **easy** instance. It must be a value between 0 and 1.
-	```opts.selvars.smallscaleflag``` by setting this flag as ```TRUE```, you can carry out a small scale experiment using a randomly selected fraction of the original data. This is useful if you have a large dataset with more than 1000 instances, and you want to explore the parameters of the model.
-	```opts.selvars.smallscale``` fraction taken from the original data on the small scale experiment.
-	```opts.selvars.fileidxflag``` by setting this flag as ```TRUE```, you can carry out a small scale experiment. This time you must provide a ```.csv``` file that contains in one column the indices of the instances to be taken. This may be useful if you want to make a more controlled experiment than just randomly selecting instances.
-	```opts.selvars.fileidx``` name of the file containing the indexes of the instances.

### Dimensionality reduction settings

The toolkit uses PBLDR as a dimensionality reduction method. Technical details about it can be found [here](https://doi.org/10.1007/s10994-017-5629-5).

-	```opts.pbldr.analytic``` determines whether the analytic (set as ```TRUE```) or the numerical (set as ```FALSE```) solution to the dimensionality reduction problem should be used. We recommend to leave this setting as ```FALSE```, due to the instability due to possible poor-conditioning of the analytical solution.
-	```opts.pbldr.ntries``` number of iterations that the numerical solution of PBLDR attempts.

###  Algorithm selection settings

The toolkit uses SVMs with radial basis kernels as algorithm selection models, through the  [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) library. To tune these models, a cross-validation experiment is carried out automatically.

-	```opts.oracle.cvgrid``` number of parameter combinations to be generated using a Latin-hypercube sampling method.
-	```opts.oracle.maxcvgrid``` maximum value for the SVM parameters. Corresponds to a power of 2.
-	```opts.oracle.mincvgrid``` minimum value for the SVM parameters. Corresponds to a power of 2.
-	```opts.oracle.cvfolds``` number of folds of the cross-validation experiment.

### Footprint construction settings

The toolkit uses Delaunay Triangulation to define the regions in the space where we statistically infer good algorithm performance. The triangulation is then pruned to remove those sections for which the evidence is poor or non-existing. Details of the algorithms employed and these settings can be found [here](https://doi.org/10.1162/evco_a_00194). We suggest not to modify the maximum and minimum sizes.

-	```opts.footprint.RHO``` minimum density required for a section of a footprint.
-	```opts.footprint.PI``` minimum purity required for a section of a footprint.
-	```opts.footprint.LOWER_PCTILE``` minimum size for the side of a triangle, as a percent of the largest distance.
-	```opts.footprint.UPPER_PCTILE``` maximum size for the side of a triangle, as a percent of the largest distance.

### Automatic data bounding and scaling

The toolkit implements simple routines to bound outliers and scale the data. **These routines are by no means perfect, and users should pre-process their data independently if preferred**. However, the automatic bounding and scaling routines should give some idea of the kind of results may be achieved. In general, we recommend that the data is transformed to become **close to normally distributed** due to the linear nature of the PBLDR optimal projection algorithm.

- ```opts.auto.preproc``` turns on (set as ```TRUE```) the automatic pre-processing.
- ```opts.bound.flag``` turns on (set as ```TRUE```) data bounding. This sub-routine calculates the median and the interquartile range ([IQR](https://en.wikipedia.org/wiki/Interquartile_range)) of each feature and performance measure, and bounds the data to the median plus or minus five times the IQR.
- ```opts.norm.flag``` turns on (set as ```TRUE```) scalling. This sub-routine scales into a positive range each feature and performance measure. Then it calculates a [box-cox transformation](https://en.wikipedia.org/wiki/Power_transform#Box%E2%80%93Cox_transformation) to stabilise the variance, and a [Z-transformation](https://en.wikipedia.org/wiki/Standard_score) to standarise the data. The result are features and performance measures that are close to normally distributed.

### Automatic feature selection

The toolkit implements simple routines to select features, given their cross-correlation and correlation to performance. Ideally, we want the smallest number of orthogonal and predictive features. **These routines are by no means perfect, and users should pre-process their data independently if preferred**.  In general, we recommend **using no more than 10 features** as input to PBLDR optimal projection algorithm, due to the numerical nature of its solution and issues in identifying meaningful linear trends.

- ```opts.auto.featsel``` turns on (set as ```TRUE```) the automatic feature selection.
- ```opts.diversity.flag``` turns on (set as ```TRUE```) the selection by diversity. This sub-routine identifies the unique values for each feature. If the fraction of unique values is below a threshold, the feature is discarded. This sub-routine works well for features that are expected to be **real numbers without much repetition between instances**. An integer/binary feature is likely to be discarded by this routine, so it should be used with caution.
- ```opts.diversity.threshold``` diversity threshold. It must be a value between 0 and 1.
- ```opts.corr.flag``` turns on (set as ```TRUE```) the selection by correlation between features and performance measures. This sub-routine calculates the [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) between the features and the performance. Then it takes its absolute value, and sorts them from largest to lowest. Then, for each algorithm, it takes the top features as defined by the threshold. It automatically bounds itself to a minimum of 3 features. We recommend using this routine, particularly if the number of features is greater than 20.
- ```opts.corr.threshold``` correlation threshold indicating the top N features per algorithm that should be selected based on correlation with performance. It should be a value higher than 1.
- ```opts.clust.flag```  turns on (set as ```TRUE```) the selection by clustering. It works by defining the [Pearson correlation coefficient](https://en.wikipedia.org/wiki/Pearson_correlation_coefficient) as a dissimilarity metric between features. Then, [k-means clustering](https://en.wikipedia.org/wiki/K-means_clustering) is used to identify groups of similar features. To select one feature per group, the algorithm first projects the subset of selected featurs into two dimensions using Principal Components Analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)) and then [Random Forests](https://en.wikipedia.org/wiki/Random_forest) to predict whether an instance is easy or not for a given algorithm. Then, the subset of features that gives the most accurate models is selected. This routine is **potentially very expensive computationally** due to the multiple layer training process. However, it is our current recommended approach to select the most relevant features. At the moment, this routine tests all possible combinations, so it should not be used if there are more than 20 features available.
- ```opts.clust.KDEFAULT``` number of default clusters. The routine assumes at least 3 clusters and no more than the number of features. Ideally it **should not** be a value larger than 10.
- ```opts.clust.SILTHRESHOLD``` the routine uses [Silohuette analysis](https://en.wikipedia.org/wiki/Silhouette_(clustering)) to automatically infer the existing number of clusters. It must be a value between 0 and 1, preferably between 0.4 and 0.8.
- ```opts.clust.NTREES``` number of threes used by the Random Forest models.
- ```opts.clust.MaxIter``` number of iterations used to converge the k-means algorithm.
- ```opts.clust.Replicates``` number of repeats carried out of the k-means algorithm.

### Output settings

These settings result in more information being stored in files or presented in the console output.

- ```opts.webproc.flag``` This flag produces the output files employed to draw the figures in MATILDA's web tools (click [here](https://matilda.unimelb.edu.au/newuser) to open an account). It is recommended to leave this setting as false.

## Contact

If you have any suggestions or ideas (e.g. for new features), or if you encounter any problems while running the code, please use the [issue tracker](https://github.com/andremun/InstanceSpace/issues) or contact us through the MATILDA's [Queries and Feedback](http://matilda.unimelb.edu.au/contact-us) page.

## Acknowledgements

Funding for the development of this code was provided by the Australian Research Council through the Australian Laureate Fellowship FL140100012.
