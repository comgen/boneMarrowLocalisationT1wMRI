This is the code that accompanies the paper - Estimating bone marrow adiposity from head MRI and identifying its genetic architecture, Kaufmann et al, eLife 2024 (https://doi.org/10.7554/eLife.101499.1)

A description of how to use this code can be found in 000_README.pdf file (https://github.com/comgen/boneMarrowLocalisationT1wMRI/blob/main/000_README.pdf)

The GWAS summary statistics are published on FUMA:
* Male and female discovery: https://fuma.ctglab.nl/browse/153909
* Male discovery: https://fuma.ctglab.nl/browse/153916
* Female discovery: https://fuma.ctglab.nl/browse/153914


Overview of the artificial neural network training, application, and validation:

a. The data used for training the network were generated by simulating 336 different head types (with an example simulation for one head type provided here), where a head type is defined by mean bone thickness, mean BM intensity and mean body fat thickness. There are 5000 datapoints per head each with 50 intensities per datapoint. OT: outer table, IT: inner table, BM: bone marrow.

b. We evaluated the performance of the neural network using two metrics: 1. the overlap between the predicted and the true BM location (defined by the simulation), and 2. the ratio between the predicted signal intensity of the BM and the true intensity of the BM.

c. Illustration of the process of extracting the intensity arrays from the MRI scan, locating the BM, and computing its intensity.

d. The algorithm was validated on a real dataset consisting of 30 individuals that were scanned 10 times each at 3-day intervals. For each scan, we computed the intensity of the outer bone (blue) and of the BM (yellow) which displayed high concordance across scans of the same individuals. The quality control measures correctly identified scans not producing reliable measures (13 of 300, 11 of which fail due to too few datapoints).

e. Application of the algorithm to a cohort of monozygotic and same-sex dizygotic twins. Correlation in BMA between twin pairs in coloured boxes.

![Overview](figure_1.png "Overview of the artificial neural network training, application, and validation")