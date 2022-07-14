function [beta, variables, pval] = mpMRI_Classifier(T) 
%% Use the two-step approach with the Lasso + Ridge regressions to 
% 1. choose the most important parameters from many mpMRI parameters using Lasso regression
% 2. find the best combination of chosen MRI parameters from step#1 using Ridge regression
%
% ------------------------------------------------------------------------------------------------------------------
% If you use this function in your work, please cite the following paper
% 
% Devan SP, Jiang X, Kang H, Luo G, Xie J, Zu Z, Stokes AM, Gore JC, McKnight CD, Kirschner AN, Xu J. Towards differentiation of brain tumor from radiation necrosis using multi-parametric MRI: Preliminary results at 4.7T using rodent models. (in review)
%
% ------------------------------------------------------------------------------------------------------------------
% INPUTS: 
%           T: a table: rows as subjects and columns as different parameters. 
%   NOTE: the last column must be the ground truth binary classifier, i.e. 0 or 1. 
% 
% ------------------------------------------------------------------------------------------------------------------
% OUTPUTS:
%       beta: the beta coefficients of chosen parameters
%       variables: the name of the chosen parameters
%
% ------------------------------------------------------------------------------------------------------------------
% EXAMPLE: 
% A T table is 
%       iAUC             ADC            T2                     APT          NOE               R1f             PSR             classifier
%     ________    _______    ________    ________    ________    _______    ________    __________
% 
%     0.067685    0.81105    0.063611    0.037705         0.05306        0.47956    0.049855               1     
%     0.051937    0.66221    0.059965    0.076146         0.089565    0.45986        0.052833              1     
%     0.066156    0.77872    0.058203    0.066225         0.090489    0.42964     0.044197                 1     
%     0.071491    0.72784    0.064424    0.088631          0.12789       0.46257    0.045446               1     
%     0.068834    0.80253    0.061734    0.073479             0.10072    0.48478    0.047072               1     
%     0.021854    0.89895    0.061944    0.033825         0.078686    0.46647      0.0417                   0     
%     0.028181        1.1       0.058994     0.02656           0.06515        0.49103    0.044089                0     
%      0.03475     1.0333    0.062767        0.028864         0.061514     0.5317    0.063902                 0     
%     0.039976     1.1246      0.0654        0.033752        0.051271       0.49942    0.035997               0     
%     0.034866    0.84665    0.056446    0.027275         0.039886       0.48131    0.042199                0     
%
% Run command:  [beta, variables, pval] = My_Lasso_Ridge_regression(T) 
% 
% Results:
% 1. The Lasso regression chooses three parameters (non-zero beta-coefficients): iAUC, ADC, and APT
% 2. The Ridge regression found eta = -1.79 + 148.77*iAUC -9.49*ADC + 61.70*APT
%
% ------------------------------------------------------------------------------------------------------------------

    %% preparation
    Data = T{:,1:(end-1)} ; 
    binary_classifier = T{:,end} ; 
    if any(~ismember(binary_classifier,[0 1]))
        error('The last column of the table T must be a binary (0 or 1) classifier') ; 
    end

    %% ------------- Step #1 ----------------
    % LASSO Logistic Regression (alpha=1)
    % meas: measured values (matrix: n_subj x n_par)
    % tissue_binary (1=tumor, 0=necrosis; row vector length n_subj)

    disp('------ Step#1: Running the Lasso regression to choose dominant parameters (with non-zero beta coefficients) ------')
    [B, FitInfo] = lassoglm(Data, binary_classifier,'binomial','Link','logit','NumLambda',100,'CV',10,'MCRep',100,'Alpha',1) ; 
    ax=lassoPlot(B,FitInfo,'PlotType','CV'); % Plot D vs. Lambda
    legend('show','Location','best') % show legend
    indx = FitInfo.Index1SE; % Choose lambda by 1 SE rule
    beta_coef = B(:,indx) ; 
    nonzeros = sum(beta_coef ~= 0) ; 
    beta0 = FitInfo.Intercept(indx) ; 
    beta = [beta0;beta_coef] ; 
    preds = glmval(beta,Data,'logit') ; 
    ind_predictors = find(beta_coef)' ;  % indices of nonzero predictors
    variables = T.Properties.VariableNames(ind_predictors) ; 
    fprintf('\t ----- Done! The Lasso regression found the following parameters... \n')
    for n=1:length(ind_predictors)
        fprintf('\t\t%d: %s\n',n, char(T.Properties.VariableNames(ind_predictors(n))))
    end
    
    %% ------------- Step #2 ----------------
    % RIDGE logistic regression (alpha -> 0)
    % Use nonzero coefficients from lasso
    disp('------ Step#2: Running the Ridge regression to find the combination of chosen parameters ------')
    [B, FitInfo] = lassoglm(Data(:,ind_predictors), binary_classifier,'binomial','Link','logit','NumLambda',100,'MCRep',100,'CV',10,'Alpha',.01);
    indx = FitInfo.Index1SE ; % Choose lambda by 1 SE rule
    beta_coef = B(:,indx) ; 
    nonzeros = sum(beta_coef ~= 0) ;  % Number of nonzero predictors
    beta0 = FitInfo.Intercept(indx) ; 
    beta = [beta0;beta_coef] ; 
    preds = glmval(beta, Data(:,ind_predictors), 'logit') ;  % Estimated P at each point
    fprintf('\t ----- Done! The Ridge regression found the following combinaiton... \n')
    str = sprintf('\t\t\t eta = %+.2f ',beta0) ; 
    for n=1:length(ind_predictors)
        str = [str, sprintf('%+.2f * %s  ', beta_coef(n), char(T.Properties.VariableNames(ind_predictors(n))))] ; 
    end
    fprintf([str, '\n']) ; 
    
    %% deviance calc adapted from the GeneralizedLinearModel class
    dev = @(mu,y)2*(y.*log((y+(y==0))./mu)+(1-y).*log((1-y+(y==1))./(1-mu))) ; 
    deviance = sum(dev(preds, binary_classifier)) ; 
    constant_model = mean(binary_classifier)*ones(size(preds)) ; 
    deviance_null = sum(dev(constant_model, binary_classifier)) ; 

    pval = 1 - chi2cdf(deviance_null-deviance,nonzeros) ; 
    fprintf('--------- The p value is %.4f ------\n', pval) ; 
    
    %% visualization
    figure('NumberTitle', 'off', 'Name', 'Fitted Ridge Regression') ; clf ; hold on ; 
    set(gcf, 'Units', 'inches', 'Position', [9 6 3.4 2.8]) ; %colormap gray

    eta = Data(:,ind_predictors)*beta_coef + beta0 ; % Linear model
    [sort_meas, i_sort] = sort(eta);

    nexp = ceil(log10(max(abs(eta)))) ; 
    x=linspace(-10^nexp,10^nexp,100) ;
    plot(eta(binary_classifier==1), binary_classifier(binary_classifier==1), 'kd', 'MarkerSize',9, 'MarkerFaceColor','r','Linewidth',1) ; 
    plot(eta(binary_classifier==0), binary_classifier(binary_classifier==0), 'bo',  'MarkerSize',9, 'MarkerFaceColor','g','Linewidth',1) ; 
    hold on
    plot(x, 1./(1+exp(-x)),'b-','LineWidth',2)
    plot([0 0],[0 1],'k--')
    set(gca,'FontSize',10, 'XLim',[-1 1]*10^nexp) ;  box on ; 
    xlabel('$\eta = \beta_0 + \sum_n\beta_nx_n$','interpreter','latex','FontSize',12)
    ylabel('Probability','FontSize',12)
    hlg = legend('9L','RN','fit') ; set(hlg, 'FontSize',11, 'color','none','edgecolor','none', 'NumColumns',1, 'location','best')
    
end
