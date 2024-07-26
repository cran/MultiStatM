## Version 2.0.0
 
* Most functions have been renamed and reorganized in their arguments and output. For example the functions \code{conv_Cum2Mom} and \code{conv_Cum2MomMulti} of version 1.2.1 have been joined in the new function \code{Cum2Mom} which now has an additional argument \code{Type=c("Univariate","Multivariate")}.

* A detailed list of the changes, connecting functions in version 1.2.1 to those of version 2.0.0 is provided in the Vignette.



## Version 1.2.1

* Minor changes for code optimization


## Version 1.2.0

* The function \code{Esti_Gram_Charlier}, to compute a Gram-Charlier approximation for a multivariate density has been introduced. 

* New functions \code{indx_Commutator_Mixing} and \code{indx_Commutator_Moment}, avoiding use of (large) matrices have been introduced.

* The codes of the following functions have been modified in order to get faster computing time: \code{Esti_MRSz_Kurt};  \code{Esti_Skew_Variance_Th}; \code{Esti_Kurt_Variance_Th}; \code{distr_Uni_EVSK_Th}; \code{distr_Uni_MomCum_Th} and \code{Hermite_N_Cov_X1_X2}.

## Version 1.1.1

* Dependence on package primes has been removed


## Version 1.1.0

* This is the first public release and distribution of the package
