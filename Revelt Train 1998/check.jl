# This code checks the input data and specifications that the user provides in mxlhb.m
# Written by Kenneth Train, July 26, 2006 (latest edits on Sept 24, 2006).
# Ported to Julia by Giuseppe Forte on Dec, 10 2020.

function check() 

    ok = 0;

    # Check for positive intergers

    if (ceil(np) != np) | (np < 1);
        println(@sprintf "\n NP must be a positive integer, but it is set to %d \n Program terminated." np);
        return
    end

    if (ceil(ncs) != ncs) | (ncs < 1);
        println(@sprintf "\n NCS must be a positive integer, but it is set to %d \n Program terminated." ncs);
        return
    end

    if (ceil(nrows) != nrows) | (nrows < 1);
        println(@sprintf "\n NROWS must be a positive integer, but it is set to %d \n Program terminated." nrows);
        return
    end

    if (ceil(ndraws) != ndraws) | (ndraws < 1);
        println(@sprintf "\n NDRAWS must be a positive integer, but it is set to %d \n Program terminated." ndraws);
        return
    end

    if (ceil(drawtype) != drawtype) | (drawtype < 1);
        println(@sprintf "\n DRAWTYPE must be a positive integer, but it is set to %d \n Program terminated." drawtype);
        return
    end

    if (ceil(seed1) != seed1) | (seed1 < 1);
        println(@sprintf "\n SEED1 must be a positive integer, but it is set to %d \n Program terminated." seed1);
        return
    end

    if (nv > 0) & sum((ceil.(idv) .!= idv) .+ (idv .< 1)) != 0;
        println(@sprintf "\n IDV must contain positive integers only, but it contains other values. \n Program terminated.");
        return
    end


    if (nf > 0) & sum((ceil.(idf) .!= idf) + (idf .< 1)) != 0;
        println(@sprintf "\n IDF must contain positive integers only, but it contains other values. \n Program terminated.");
        return
    end

    # Checking XMAT #
    if size(xmat, 1) != nrows
        println(@sprintf "\n XMAT has %d rows, but it should have %d rows. \n Program terminated." size(xmat, 1) nrows);
        return
    end

    if sum(xmat[:, 1] .> np) != 0
        println(@sprintf "\n The first column of XMAT has a value greater than NP = %d. \n Program terminated." np);
        return
    end

    if sum(xmat[:, 1] .< 1) != 0
        println(@sprintf "\n The first column of XMAT has a value smaller than 1. \n Program terminated.");
        return
    end

    k = (xmat[2:nrows, 1] .!= xmat[1:nrows - 1, 1]) .& (xmat[2:nrows, 1] .!= (xmat[1:nrows - 1, 1] .+ 1.0));
    if sum(k) != 0
        println("\n The first column of XMAT does not ascend from 1 to NP.\n Program terminated.");
        return
    end

    if sum(xmat[:, 2] .> ncs) != 0
        println(@sprintf "\n The second column of XMAT has a value greater than NCS = %d. \n Program terminated." ncs);
        return
    end

    if sum(xmat[:, 2] .< 1) != 0
        println("\n The second column of XMAT has a value less than 1. \n Program terminated.");
        return
    end

    k = (xmat[2:nrows, 2] .!= xmat[1:nrows - 1, 2]) .* (xmat[2:nrows, 2] .!= (xmat[1:nrows - 1, 2] .+ 1.0));
    if sum(k) != 0
        println("\n The second column of XMAT does not ascend from 1 to NCS. \n Program terminated.")
        return
    end

    if sum((xmat[:, 3] .!= 0.0) .& (xmat[:, 3] .!= 1.0)) != 0
        println("\n The third column of XMAT has a value other than 1 or 0. \n Program terminated.")
        return
    end

    for s = 1:ncs
        k = xmat[:, 2] .== s;
        if sum(xmat[k, 3]) > 1
            println(@sprintf "\n The third column of XMAT indicates more than one chosen alternative for choice situation %d. \n Program terminated." s);
            return
        end
        if sum(xmat[k, 3]) < 1
            println(@sprintf "\n The third column of XMAT indicates that no alternative was chosen for choice situation %d. \n Program terminated." s);
            return
        end 
    end

    if sum(isnan.(xmat)) != 0
        println("\n XMAT contains missing data. \n Program terminated.");
        return
    end;

    if sum(isinf.(xmat)) != 0
        println("\n XMAT contains an infinite value. \n Program terminated.");
        return
    end;

    if (nv > 0) & (size(idv, 2) != 2);
        println(@sprintf "\n IDV must have 2 columns and yet it is set to have %d. \n Program terminated." size(idv, 2));
        return
    end;

    if (nv > 0) & (sum(idv[:, 1] .> size(xmat, 2)) != 0);
        println("\n IDV identifies a variable that is outside XMAT. \n The first column of IDV is");
        idv[:, 1]'
        println(@sprintf "\n when each element of this row must be no greater than %d, which is the number of columns in XMAT. \n Program terminated." size(XMAT,2))
        return
    end;

    if (nv > 0) & (sum(idv[:, 1] .≤ 3) != 0);
        println("\n Each element in the first column of IDV must exceed 3, since the first three variables in XMAT cannot be explanatory variables.");
        println("\n But the first column of IDV is");
        idv[:, 1]'
        println("\n which has an element below 3. \n Program terminated.");
        return
    end;

    if (nv > 0) & (sum((idv[:, 2] .< 0) .+ (idv[:, 2] .> 6)) != 0);
        println("\n The second column of IDV must be integers 1-6 identifying the distributions. But the second column of IDV is specified as");
        idv[:, 2]'
        println("\n which contains a number other than 1-6. \n Program terminated.");
        return
    end;

    if (nv > 0) & (size(names_idv, 2) != 1);
        println(@sprintf "\n NAMES must have 1 column and yet it is set to have %d. Be sure to separate names by semicolons. \n Program terminated." size(names_idv, 2));
        return
    end;

    if (nv > 0) & (size(idv, 1) != size(names_idv, 1));
        println(@sprintf "\n IDV and NAMES must have the same length but IDV has length %d, while NAMES has length %d. \n Program terminated." size(idv, 1) size(names_idv, 1));
        return
    end; 

    if (nv > 0) & (size(B, 2) != 1);
        println(@sprintf "\n B must have 1 column and yet it is set to have %d. Be sure to separate values by semicolons. \n Program terminated." size(B, 2));
        return
    end;
    
    if (nv > 0) & (size(B, 1) != size(idv, 1));
        println(@sprintf "\n B must have the same length as IDV but instead has length %d. \n Program terminated." size(B, 1));
        return
    end; 

    if (nv > 0) & (size(W, 2) != 1);
        println(@sprintf "\n W must have 1 column and yet it is set to have %d. Be sure to separate values by semicolons. \n Program terminated." size(W, 2));
        return
    end;
    
    if (nv > 0) & (size(W, 1) != size(idv, 1));
        println(@sprintf "\n W must have the same length as IDV but instead has length %d. \n Program terminated." size(W, 1));
        return
    end; 

    if (nf > 0) & (size(idf, 2) != 1);
        println(@sprintf "\n IDF must have 1 column and yet it is set to have %d. Be sure to separate elements by semicolons. \n Program terminated." size(idf, 2));
        return
    end;

    if (nf > 0) & (sum(idf .> size(xmat, 2)) != 0);
        println("\n IDF identifies a variable that is outside XMAT. IDF is");
        idf'
        println(@sprintf "\n where each element must be no greater than %d which is the number of columns in XMAT. \n Program terminated." size(xmat, 2));
        return
    end;

    if (nf > 0) & (sum(idf .≤ 3) != 0);
        println("\n Each element of IDF must exceed 3 since the first three variables of XMAT cannot be explanatory variables. But IDV is");
        idv[:, 1]'
        println("\n which contains an element below 3. \n Program terminated.")
        return
    end;

    if (nf > 0) & (size(names_idf, 2) != 1);
        println(@sprintf "\n NAMES_IDF must have 1 column and yet it is set to have %d. Be sure to separate names by semicolons. \n Program terminated." size(names_idf, 2));
        return
    end;

    if (nf > 0) & (size(idf, 1) != size(names_idf, 1));
        println(@sprintf "\n IDF and NAMESF must have the same length but IDF has length %d, while NAMESF has length %d. \n Program terminated." size(idf, 1) size(names_idf, 1));
        return
    end; 

    if (nf > 0) & (size(F, 2) != 1);
        println(@sprintf "\n F must have 1 column and yet it is set to have %d. Be sure to separate values by semicolons. \n Program terminated." size(F, 2));
        return
    end;
    
    if (nf > 0) & (size(F, 1) != size(idf, 1));
        println(@sprintf "\n F must have the same length as IDF but instead has length %d. \n Program terminated." size(F, 1));
        return
    end; 

    if sum(drawtype .< 0 | drawtype .> 4) != 0;
        println(@sprintf "\n DRAWTYPE must be an integer 1-4 identifying the type of draws. But DRAWTYPE is set to %d. \n Program terminated." drawtype);
        return
    end;

    if wantwgt != 0 & wantwgt != 1
        println(@sprintf "\n WANTWGT must be 0 or 1, but it is set to %d. \n Program terminated." wantwgt);
        return
    end

    if (wantwgt == 1) & (size(idwgt, 1) != 1)
        println("\n When WANTWGT==1, as you have set it, then IDWGT must be an scalar that identifies a variable in XMAT for the weights.");
        println(@sprintf "\n But IDWGT is set to %d. \n Program terminated." idwgt);
    return
    end

    if (wantwgt == 1) & (idwgt > size(xmat, 2))
        println("\n When WANTWGT==1, as you have set it, then IDWGT must identify a variable in XMAT for the weights.");
        println(@sprintf "\n But IDWGT is set to %d when XMAT has only %d columns. \n Program terminated." idwgt size(xmat, 2));
        return
    end

    if (wantwgt == 1) & (idwgt < 1 | ceil(idwgt) != idwgt)
        println(@sprintf "\n IDWGT must be a positive integer indentifying a variable in XMAT but it is set to %d. \n Program terminated." idwgt);
        return
    end

    if wantwgt == 1
        cp = xmat[:, 1];
        for r = 1:np
            if sum(xmat[cp .== r, idwgt] .!= mean(xmat[cp .== r, idwgt])) > 0
                println(@sprintf "\n Variable identified by IDWGT is %d. This weight variable must be the same for all rows of data for each person." idwgt);
                println(@sprintf "\n However, it is not the same for all rows for person %d, and maybe for people after that person (not checked). \n Program executed." r);
                return
            end
        end
    end

    ok = 1;
end 