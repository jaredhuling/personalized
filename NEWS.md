# personalized 0.2.2


# personalized 0.2.1

* Fixes default argument bug for r-oldrel + windows

# personalized 0.2.0

* Simplified plot labeling
* Added clarifications to documentation
* Added customized loss function option
* Added options for count outcomes via Poisson negative log-likelihood as the loss
* Added treatment effect calculation based on estimated benefit scores
* Clarified/improved printing
* Improved numerical stability of weighted.ksvm

# personalized 0.1.5

* Added printing of subgroup_validated results for quantile cutoffs via the which.quant argument

# personalized 0.1.4

* Added plots of means within treatment groups as the benefit scores are varied
* Added quantile and median cutpoints as options
* Fixes subgroup effect calculation to account for weights

# personalized 0.1.3

* Added NSW Study dataset

# personalized 0.1.2

* Added requirement for latest version of glmnet. old versions throw error when efficiency augmentation used

# personalized 0.1.1

* Fixed minor bugs regarding multiple treatment options, match.id

* Added OWL-type losses: logistic and hinge surrogates (multiple treatment available for logistic surrogate)

* Added outcome flipping OWL-type losses: logistic and hinge surrogates

* Added augmentation option for non-continuous outcomes via offset

# personalized 0.1.0

* Added estimation functionality for multiple treatments

* Updated all plot and summary type functions to properly handle results for multiple treatments

* Aaron added plotting option for validation objects that allows the user to inspect the distribution of variable selections over the bootstrap or training/test resampling iterations

* Added more options to printing of validation objects

* Updated check.overlap() to handle multiple treatments

* Aaron added match.id argument to allow proper analysis of matched case-control datasets
