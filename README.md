# MATILDA

Input parameters:

X:    A cell matrix of features. The first row corresponds to the feature names, while the remainder should have the data values.

Y:    A cell matrix of algorithm performance. The first row corresponds to the algorithm names, while the remainder should have the data values.

Ybin: A cell matrix with a binary measure of algorithm performance, where '1' is good performance and '0' is bad. The first row corresponds to the algorithm names.

opts: An structure with all the options. These can be found in the 'example.m' file.

Output paramters:

out: An structure with all the outputs. The most important results are those in out.footprint.performance which contains the table with areas and densities, and out.pbldr.A that contains the projection matrix. The location of the instances in the space is defined by out.pbldr.Z
