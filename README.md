# Estimators for Discrete Choice models, collected as I go. 

*Please direct comments/questions about the code to me, not to the authors of the papers I am citing.*

### Revelt Train (1998)
This is a straightforward port of Kenneth Train's [Mixed Logit Maximum Simulated Likelihood](https://eml.berkeley.edu/Software/abstracts/train1006mxlmsl.html) code
as detailed in Revelt Train (1998; publisher version [here](https://www.mitpressjournals.org/doi/10.1162/003465398557735), working paper [here](https://eml.berkeley.edu/wp/train0797b.pdf)). 

Modifications of note: 
- Removed the possibility to use stored quasirandom draws;
- Removed the possibility to store in memory only a subset of random draws at a time;
- Added the possibility of not using analytical gradient in likelihood (though it is advisable to do so);
- The likelihood/gradient function is modified to run faster; more comments add context to commands.

To run: parameters are defined and data is imported by `mxlmsl.jl`. Note the `data.txt` called in the file is available from [Train's original upload](https://eml.berkeley.edu/Software/source_code/train_mxlmsl_06.zip). `loglik.jl` and `llgrad2.jl` are included but the script only calls `llg.jl`, which supersedes the aforementioned files.

### Train (2016)
This is a straightforward port of Kenneth Train's [Mixed logit with a flexible mixing distribution](https://eml.berkeley.edu/~train/flexsplash.html) code
as detailed in Train (2016; publisher version [here](https://eml.berkeley.edu/~train/flexible.pdf)). This paper generalises ideas found in Section 6 of Train (2008; publisher version [here](https://eml.berkeley.edu/~train/EMtrain.pdf)).

Modifications of note: 
- Removed the possibility to use GPU computing (could be added, I am not acquainted with Julia possibilities);
- Added the possibility to estimate the model in preference space;
- Added the possibility to not subset from S to S_r (see Train 2016, p. 42).

To run: parameters are defined and data is imported by `FlexibleWtp.jl`. Note the `videodata100.mat` called in the file is available from [Train's original upload](http://eml.berkeley.edu/~train/flexcodes.zip). Estimation is performed by `doestimation.jl`, and bootstrap is available via `boot.jl`. 
