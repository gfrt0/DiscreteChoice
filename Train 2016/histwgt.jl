# Function to create histogram bins of distribution approximations
# Written in Matlab by Kenneth Train, first version Oct 21, 2015 (correction to file CreateZ.m made on Feb 10, 2020).
# Ported to Julia by Giuseppe Forte Dec, 22, 2020.

function histwgt(vv, ww, startpt, endpt, nbins) 

    #Inputs: 
    #vv: column vector of values 
    #ww: column vestor of weights 
    #startpt - lowest value to use for histogram 
    #endpt - highest value to use for histogram.
    #nbins - number of bins 
    #Histogram has nbins evenly spaced bins between startpt and endpt. 

    #Outputs: 
    #histwg: nbinsx1 vector: share of weights in each bin 
    #midbin: nbinsx1: midpoint value for each bin

    δ = (endpt - startpt) / nbins; 
    subs = ceil.((vv .- startpt) / δ); 
    subs[subs .== 0] .= 1;
    midbin = ((startpt + δ / 2):δ:endpt)';
    histwg = [sum(ww[subs .== j]) for j in 1:nbins]

    return histwg, midbin
end


 


 
