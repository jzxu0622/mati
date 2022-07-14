%% Eg05_mpMRI_Classifier
% This example demonstrate how to use the function mpMRI_Classifier.m with Lasso and Ridge regressions to 
% 1. choose the most important parameters from many mpMRI parameters using Lasso regression
% 2. find the best combination of chosen MRI parameters from step#1 using Ridge regression
%
% ------------------------------------------------------------------------------------------------------------------
% If you use this function in your work, please cite the following paper
% 
% Devan SP, Jiang X, Kang H, Luo G, Xie J, Zu Z, Stokes AM, Gore JC, McKnight CD, Kirschner AN, Xu J. Towards differentiation of brain tumor from radiation necrosis using multi-parametric MRI: Preliminary results at 4.7T using rodent models. (in review)
%
% ------------------------------------------------------------------------------------------------------------------

%% load condensed data (tumor vs RN)
load Eg05_mpMRI_9L_RN.mat

[beta, variables, pval] = mpMRI_Classifier(T) ; 

disp('End of the program Eg05_mpMRI_Classifier.m') 
disp('-------- Update the figure labels/legends/axes/... according to your needs ... ------') ; 


