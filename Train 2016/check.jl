# This code checks the input data and specifications that the user provides in FlexibleWtp.m

# Check for positive intergers

function check()

   check_ok = false;

   if ceil(NP) != NP | NP < 1;
      println(@sprintf "NP must be a positive integer, but it is set to %d. \nProgram terminated.\n" NP);
      return
   end

   if ceil(NCS) != NCS | NCS < 1;
      println(@sprintf "NCS must be a positive integer, but it is set to %d. \nProgram terminated.\n" NCS);
      return
   end

   if ceil(NROWS) != NROWS | NROWS < 1;
      println(@sprintf "NROWS must be a positive integer, but it is set to %d. \nProgram terminated.\n" NROWS);
      return
   end

   # Checking XMAT #
   if (size(XMAT,1) != NROWS)
      println(@sprintf "XMAT has %d rows, but it should have NROWS = %d rows. \nProgram terminated.\n" size(XMAT, 1) NROWS);
      return
   end

   if sum(XMAT[:, 1] .> NP) != 0
      println(@sprintf "The first column of XMAT has a value greater than NP = %d. \nProgram terminated.\n" NP);
      return
   end

   if sum(XMAT[:,1] .< 1) != 0
      println("The first column of XMAT has a value less than 1. \nProgram terminated.\n");
      return
   end

   k = (XMAT[2:NROWS, 1] .!= XMAT[1:NROWS - 1, 1]) .& (XMAT[2:NROWS, 1] .!= (XMAT[1:NROWS - 1, 1] .+ 1));
   if sum(k) != 0
      println("The first column of XMAT does not ascend from 1 to NP. \nProgram terminated.\n");
      return
   end

   if sum(XMAT[:, 2] .> NCS) != 0
      println(@sprintf "The second column of XMAT has a value greater than NCS = %d. \nProgram terminated.\n" NCS);
   end

   if sum(XMAT[:, 2] .< 1) != 0
      println("The second column of XMAT has a value less than 1. \nProgram terminated.\n");
      return
   end

   k = (XMAT[2:NROWS, 2] .!= XMAT[1:NROWS - 1, 2]) .& (XMAT[2:NROWS, 2] .!= (XMAT[1:NROWS - 1, 2] .+ 1));
   if sum(k) != 0
      println("The second column of XMAT does not ascend from 1 to NCS. \nProgram terminated.\n");
      return
   end


   if sum(XMAT[:, 3] .!= 0 .+ XMAT[:, 3] .!= 1) != 0
      println("The third column of XMAT has a value other than 1 or 0. \nProgram terminated.\n");
      return
   end

   for s = 1:NCS
      kk = (XMAT[:, 2] .== s);
      if sum(XMAT[kk, 3]) > 1
         println(@sprintf "The third column of XMAT indicates more than one chosen alternative for choice situation %d. \nProgram terminated.\n" s);
         return
      end
      if sum(XMAT[kk, 3]) < 1
         println(@sprintf "The third column of XMAT indicates no chosen alternative for choice situation %d. \nProgram terminated.\n" s);
         return
      end 
   end

   if sum(isnan.(XMAT)) != 0
      println("XMAT contains missing data. \nProgram terminated.\n");
      return
   end;

   if sum(isinf.(XMAT)) != 0
      println("XMAT contains an infinite value. \nProgram terminated.\n");
      return
   end;

   if ceil(IDPRICE) != IDPRICE | IDPRICE < 4;
      println("IDPRICE must be a positive integer >3, because the first three variables in XMAT cannot be explanatory variables.");
      println(@sprintf "But IDPRICE is set to %d. \nProgram terminated.\n" IDPRICE);
      return
   end

   if IDPRICE > size(XMAT, 2);
      println("IDPRICE is set to a number greater than the columns in XMAT. \nProgram terminated.\n");
      return
   end

   if sum(IDV[:, 1] .> size(XMAT,2)) != 0;
      println("IDV identifies a variable that is outside XMAT.");
      println("IDV is");
      IDV[:, 1]
      println(@sprintf "when each element of this vector must be no greater than %d, the number of columns in XMAT. \nProgram terminated.\n" size(XMAT,2));
      return
   end;

   if sum(IDV[:, 1] .<= 3) != 0;
      println("Each element of IDV must exceed 3 since the first three variables in XMAT cannot be explanatory variables.");
      println("But IDV is");
      IDV[:, 1]
      println("which has an element below 3. \nProgram terminated.\n")
      return
   end;

   if size(NAMES, 2) != 1;
      println(@sprintf "NAMES must have 1 column and yet it is set to have %d." size(NAMES, 2));
      println("Be sure to separate names by semicolons. \nProgram terminated.\n");
      return
   end;

   if size(NAMES, 1) != NV;
      println(@sprintf "NAMES must equal length(IDV) + 1 but it has length %d, while length(IDV) + 1 = %d. \nProgram terminated.\n" size(NAMES, 1) (length(IDV)+1));
      return
   end; 

   if size(P_Range, 2) != 2;
      println(@sprintf "P_Range must have 2 columns and yet it is set to have %d. \nProgram terminated.\n" size(P_Range, 2));
      return
   end;

   if size(P_Range, 1) != 1;
      println(@sprintf "P_Range must have 1 row and yet it is set to have %d. \nProgram terminated.\n" size(P_Range, 1));
      return
   end;

   if P_Range[1, 1] >= P_Range[1, 2];
      println("The second element of P_Range must exceed the first. \nProgram terminated.\n");
      return
   end;

   if size(X_Range, 2) != 2;
      println(@sprintf "X_Range must have 2 columns and yet it is set to have %d. \nProgram terminated.\n" size(X_Range, 2));
      return
   end;

   if size(X_Range, 1) != size(IDV, 1);
      println("X_Range must have the same length as IDV. \nProgram terminated.\n");
      return
   end;

   if sum(X_Range[:, 1] .>= X_Range[:, 2]) >0;
      println("The second column of X_Range must exceed the first column in all rows. \nProgram terminated.\n");
      return
   end;


   if ceil(NGridPts) != NGridPts | NGridPts < 1;
      println(@sprintf "NGridPts must be a positive integer, but it is set to %d. \nProgram terminated.\n" NGridPts);
      return
   end

   if ceil(NDRAWS) != NDRAWS | NDRAWS < 1;
      println(@sprintf "NDRAWS must be a positive integer, but it is set to %d. \nProgram terminated.\n" NDRAWS);
      return
   end

   if ceil(ThisSeed) != ThisSeed | ThisSeed < 1;
      println(@sprintf "ThisSeed must be a positive integer, but it is set to %d. \nProgram terminated.\n" ThisSeed);
      return
   end

   if ZTYPE != 1 & ZTYPE != 2 & ZTYPE != 3 ;
      println(@sprintf "ZTYPE must be 1,2, or 3, but it is set to %d. \nProgram terminated.\n" ZTYPE);
      return
   end

   if ZTYPE==1
      if ceil(PolyOrder) != PolyOrder | PolyOrder < 1;
         println(@sprintf "PolyOrder must be a positive integer, but it is set to %d. \nProgram terminated.\n" PolyOrder);
         return
      end
   end

   if ZTYPE==2
      if ceil(NLevels) != NLevels | NLevels < 1;
         println(@sprintf "NLevels must be a positive integer, but it is set to %d. \nProgram terminated.\n" NLevels);
         return
      end
   end

   if ZTYPE==3
      if ceil(NKnots) != NKnots | NKnots < 1;
         println(@sprintf "NKnots must be a positive integer, but it is set to %d." NKnots);
         return
      end
   end

   if CrossCorr ∉ [0, 1];
      println(@sprintf "CrossCorr must be 0 or 1, but it is set to %d. \nProgram terminated.\n" CrossCorr);
      return
   end


   if size(StartB, 1) != NZ;
      println(@sprintf "StartZ must contain NZ elements in a column vector, but it %d rows. \nProgram terminated.\n" size(StartB, 1));
      return
   end

   if size(StartB,2) != 1;
      println(@sprintf "StartB must have 1 column and yet it is set to have %d. " size(StartB, 2));
      println("Be sure to separate values by semicolons. \nProgram terminated.\n");
      return
   end;

   if ceil(NBins) != NBins | NBins < 1;
      println(@sprintf "NBins must be a positive integer, but it is set to %d. \nProgram terminated.\n" NBins);
      return
   end

   # if WantHessian ∉ [0, 1];
   #    println(@sprintf "WantHessian must be 0 or 1, but it is set to %d. \nProgram terminated.\n" WantHessian);
   #    return
   # end

   if WantBoot ∉ [0, 1];
      println(@sprintf "WantBoot must be 0 or 1, but it is set to %d. \nProgram terminated.\n" WantBoot);
      return
   end

   if WantBoot == 1
      if ceil(NReps) != NReps| NReps < 1;
         println(@sprintf "NReps must be a positive integer, but it is set to %d. \nProgram terminated.\n" NReps);
         return
      end
   end

   # if YesGPU ∉ [0, 1];
   #    println(@sprintf "YesGPU must be 0 or 1, but it is set to %d. \nProgram terminated.\n" YesGPU);
   #    return
   # end

   if subsetS ∉ [0, 1];
      println(@sprintf "subsetS must be 0 or 1, but it is set to %d. \nProgram terminated.\n" subsetS);
      return
   end

   if subsetS == 0 & NDRAWS < NGridPts;
      println(@sprintf "If subsetS = 0, you must specify NDRAWS ≥ NGridPts (preferably a direct multiple, or this will be enforced). \nProgram terminated.\n" subsetS);
      return
   end

   if ceil(MAXITERS) != MAXITERS| MAXITERS < 1;
      println(@sprintf "MAXITERS must be a positive integer, but it is set to %d. \nProgram terminated.\n" MAXITERS);
      return
   end

   check_ok = true;
end 
  


