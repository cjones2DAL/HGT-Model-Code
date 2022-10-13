# HGT-Model-Code
Contains MATLAB scripts for a horizontal gene transfer model.

Figure01_Code.m - This code runs the deterministic model to generate the data used to produce Fig. 1 in the in article. 

The relevant parameters are:

tmax  = the number of model iterations or "mappings" to be run (1 x 10^4)

delta = the probability that the environment will shift to the neutral state per model iteration (0.01); 

beta1 and beta2 = two values for the level of opportunity to enter naive microbial populations (0.06 and 0.08); 

Nmax = the maximum size of the metapopulation of transferable genes (1 x 10^4).

Figure02_Code.m - This code runs the stochastic model to generate the data used to produce Fig. 2 in the in article. 
The relevant parameters are the same as those used to produce Fig. 1, with the exception of tmax, which was set to 2 x 10^4.

Figure03_Code.m - This code runs the stochastic model using different value for beta in {0.02, 0.04, 0.08, 0.10, 0.12}. 
The model is run nTrials = 50 times for each value of beta. The resuling data is then used to produce Fig. 3 in the article.
RunDeterministicModel.m and RunStochasticModel.m are MATLAB functions that are called in one or more of the above three scripts.
