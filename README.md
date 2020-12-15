# Estimators for Discrete Choice models, collected as I go. 

### Revelt Train (1998)
This is a straightforward port of Kenneth Train's [Mixed Logit Maximum Simulated Likelihood](https://eml.berkeley.edu/Software/abstracts/train1006mxlmsl.html) code
as detailed in Revelt Train (1998; publisher version [here](https://www.mitpressjournals.org/doi/10.1162/003465398557735), working paper [here](https://eml.berkeley.edu/wp/train0797b.pdf)). 

Modifications of note: 
- Removed the possibility to use stored quasirandom draws;
- Removed the possibility to store in memory only a subset of random draws at a time;
- Added the possibility of not using analytical gradient in likelihood (though it is advisable to do so);
- The likelihood/gradient function is modified to run faster; more comments add context to commands.

To run: parameters are defined and data is imported by `mxlmsl.jl`. Note the `data.txt` called in the file is available [in Train's original upload](https://eml.berkeley.edu/Software/source_code/train_mxlmsl_06.zip).
