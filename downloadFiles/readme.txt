NKFD software
Version 0.1		Wednesday 12 Dec 2007 at 05:12

This software release resulted from a request for access to the software from the ICML 2001 paper. Currently not all the experiments are recreateable because I no longer have access to all the data that Bernhard and I used. 


MATLAB Files
------------

Matlab files associated with the toolbox are:

nkfdInvGetNormAxesPoint.m: Take a point on a plot and return a point within the figure.
nlpost.m: Computes the posteriors for a noiseylabel model
demNkfd1.m: Demonstration of noise model with Gaussian distributions. First demo in ICML paper.
nlem.m: Noisy label EM algorithm for full covariance Gaussian model.
nkfdPost.m: Computes the class posterior probabilities of a KFD Model.
nkfdToolboxes.m: Load in the relevant toolboxes for noisy kernel fishers discriminant.
demNkfd2.m: Demonstration of noise model with Gaussian distributions. 
nkfdEm.m: Noisy Kernel Fisher Discriminant EM algorithm.
demNkfdFourNine.m: Compare classification of USPS 4 vs 9 for a range of noise values.
nkfdClassVisualise.m: Sets up visualisation of decision boundary for NKFD model.
nkfdActiv.m: Computes the activations of a KFD Model.
