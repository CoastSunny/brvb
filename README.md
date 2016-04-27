# brvb

 Backward Renormalization Priors and the Cortical Source Localization Problem with EEG or MEG

Matlab code implementing EEG source localization from scalp activity, using Variational Bayes:

http://www.scholarpedia.org/article/Source_localization#Sparse_Bayesian_learning_.28SBL.29_and_automatic_relevance_determination_.28ARD.29)

in multiple scales, using the variance after conversion in each scale as prior information to start the next better informed. 

The goal is to localize the sources in a fine grained scale (~2000 to ~10000 dipoles) as usual, but starting from a better starting point compared to constant prior variance accross all dipoles.

For mor information:

http://arxiv.org/abs/1502.03481
