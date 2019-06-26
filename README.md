Instance Space Analysis: A toolkit for the assesment of algorithm power
===

Instance Space Analysis is a methodology for the assesment of the strengths and weaknesses of an algorithm, and thus, an approach to objectivelly compare algorithm power. At its core its the modelling of the relationship between structural properties of an instance to the performance of a group of algorithms. Instance Space Analysis allows the construction of **footprints** for each algorithm, defined as regions in the space where we statistically infer good performance. Other insights that can be gathered from the Instance Space Analysis are:

* Visualise the distribution and diversity of existing benchmark and real-world instances.
* Assess the adequacy of the features used to characterise an instance
* Partition the instance space into recommended regions for automated algorithm selection
* Distinguish areas of the instance space where it may be useful to generate additional instances to gain further insights.

The unique advantage of visualizing algorithm performance in the instance space, rather than as a small set of summary statistics averaged across a large collection of instances, is the nuanced analysis that becomes possible to examine interesting variations in performance that may be hidden by tables of summary statistics.

This repository provides a set of MATLAB tools to carry out a complete Instance Space Analysis in an automated pipeline. It is also the  computational engine that powers the Melbourne Algorithm Test Instance Library with Data Analytics ([MATILDA](http://matilda.unimelb.edu.au/)) web tools. For further information on the Instance Space Analysis methodology can be found [here](http://matilda.unimelb.edu.au/our-methodology).


Installation Instructions
---

The only requirement for the software to run is to have a current version of MATLAB. It has been tested and known to work properly in version r2018b. Earlier versions may fail to support the processing of ```.json``` files. To make use of all functionalities, the [LIBSVM](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) library must be installed. This repository contains the Windows binaries of LIBSVM, but other platforms require compilation from source code. We suggest to add the software folder to the MATLAB path.


Working with the code
---

The main interfase is the script ```example.m``` which provides the path for the ```metadata.csv``` file, and constructs the ```options.json``` file. The path provided will also be the location of all the software outputs, such as images (in ```.png``` format), tables (in ```.csv``` format) and raw intermediate data (in ```.mat``` format). The ```metadata.csv``` file should contain a table where each row corresponds to a problem instance, and each column must strictly follow the naming convention mentioned below:

* **instances** instance identifier - We expect instance identifier to be of type "String". This column is mandatory.
* **source** instance source - This column is optional
* **feature_name** The keyword "feature_" concatenated with feature name. For instance, if feature name is "density", header name should be mentioned as "feature_density". If name consists of more than one word, each word should be separated by "_" (spaces are not allowed). There **must be more than two features** for the software to work. We expect the features to be of the type "Double".
* **algo_name** The keyword "algo_" concatenated with algorithm name. For instance, if algorithm name is "Greedy", column header should be "algo_greedy". If name consists of more than one word, each word should be separated by "_" (spaces are not allowed). You can add the performance of more than one algorithms in the same csv. We expect the algorithm performance to be of the type "Double".

Moreover, emtpty cells, NaN or null values are not allowed. We expect you to handle missing values in your data before processing. You may use [this file](https://matilda.unimelb.edu.au/matildadata/graph_coloring_problem/metadata.csv) as reference.


Options
---

The script ```example.m``` contains code that prepares all the settings used by the code. Broadly, there are settings required for the analysis itself, and settings for the pre-processing of the data. For the former these are:

* ```opts.perf.MaxMin``` determines whether the algorithm performance values provided are **efficiency** measures (set as ```TRUE```), that should be maximised, or **cost** measures that should be minimised (set as ```FALSE```).
* ```opts.perf.AbsPerf``` determines whether good performance is defined absolutelly, e.g., missclassification error is lower than a 20%, (set as ```TRUE```), or if it is defined relativelly to the best performing algorithm, e.g., missclassification error is within at least 5% of the best algorithm, (set as ```FALSE```).
* ```opts.perf.epsilon``` corresponds to the threshold used to calculate good performance. It must be of the type "Double".
* ```opts.general.betaThreshold``` corresponds to the fraction of algorithms in the portfolio that must have good performance in the instance, for it to be considered an **easy** instance. It should be a value between 0 and 1.
* ```opts.pbldr.analytic``` determines whether the analytic (set as ```TRUE```) or the numerical solution (set as ```FALSE```) solution to the dimensionality reduction problem should be used. Technical details about our custom dimensionality reduction algorithm PBLDR can be found [here](https://doi.org/10.1007/s10994-017-5629-5). We recommend to leave this setting as ```FALSE```, due to the instability of the analyitical solution.
* ```opts.pbldr.ntries``` number of iterations that the numerical solution of PBLDR is attempted.
* ```opts.oracle.cvgrid``` number of parameter combinations used to tune the algorithm selection models. These are SVM models with radial basis kernels and two parameters. The parameters are generated using a Latin-hypercube sampling method.
* ```opts.oracle.maxcvgrid``` maximum value for the SVM parameters. Corresponds to a power of 2.
* ```opts.oracle.mincvgrid``` minimum value for the SVM parameters. Corresponds to a power of 2.
* ```opts.oracle.cvfolds``` number of folds of the cross-validation experiment used to tune the parameters of the algorithm selection models.
* ```opts.footprint.RHO``` minimum density required for a section of a footprint. Footprint analysis constructs a poligon using Delaunay Triangulation. To determine whether each triangle belongs or not to the footprint, 
* ```opts.footprint.PI```
* ```opts.footprint.LOWER_PCTILE```
* ```opts.footprint.UPPER_PCTILE```
* ```opts.selvars.smallscaleflag``` by setting this flag as ```TRUE```, you can carry out a small scale experiment using a randomly selected fraction of the original data. This is useful if you have a large dataset with more than 1000 instances, and you want to set up the parameters of the model.
* ```opts.selvars.smallscale``` fraction taken from the original data on the small scale experiment.
* ```opts.selvars.fileidxflag``` by setting this flag as ```TRUE```, you can carry out a small scale experiment. This time you must provide a ```.csv``` file that contains in one column the indexes of the instances to be taken. This may be useful if you want to make a more controlled experiment than just randomly selecting instaces.
* ```opts.selvars.fileidx``` name of the file containing the indexes of the instances.

The software follows a pipeline that attempts to prepare the data for the dimensionality reduction and footprint calculation stages. 

* ```opts.auto.preproc```
* ```opts.bound.flag```
* ```opts.norm.flag```
* ```opts.auto.featsel```
* ```opts.diversity.flag```
* ```opts.diversity.threshold```
* ```opts.corr.flag```
* ```opts.corr.threshold```
* ```opts.clust.flag```
* ```opts.clust.KDEFAULT```
* ```opts.clust.SILTHRESHOLD```
* ```opts.clust.NTREES```
* ```opts.clust.MaxIter```
* ```opts.clust.Replicates```

Some auxiliary options are

* ```opts.webproc.flag```


Contact
-------

If you have any suggestions or ideas (e.g. for new features), or if you encounter any problems while running the code, please use the issue tracker or contact us through the MATILDA's [Queries and Feedback](http://matilda.unimelb.edu.au/contact-us) page.



