function [vdist, resnorm] = FitFcn(this, signal, noiseSigma)
%% MRICytometry.FitFcn() fits data to MRI-cytometry model for cell size distribution
% INPUTS: 
%       this:                   MRICytometry object
%       signal:               measured signal. should be a Nacq X 1 vector
%       noiseSigma:    noise sigma
% OUTPUTS: 
% parmap is a Nx x Ny x Nz cell array. Each cell contains
%       parmap{1}: dist of intra dvw  of cell-volume-weighted (cell size)
%       parmap{2}: dist of intra Din of cell-volume-weighted
%       parmap{3}: dist of extra Dex
%       parmap{4}: dist of extra betaex
%       parmap{5}: dist of intra d non-volume-weighted (cell size)
%       parmap{6}: dist of intra Din of non-volume-weighted
%       parmap{7}: dist of intra dCal non-volume-weighted (cell size), calculated from parmap{1} to show lower-pass filter effects
%       parmap{8}: dist of intra Din non-volume-weighted, calculated from parmap{1}
%       parmap{9}: TBD
%       parmap{10}: TBD
%       parmap{11}: matrixInvw (cell-volume-weighted) obtained in step#1
%       parmap{12}: matrixExvw (cell-volume-weighted) obtained in step#1
%       parmap{13}: matrixIn (non-cell-volume-weighted) obtained in step#2
%       parmap{14}: matrixInCal (none-cell-volume-weighted) calculated directly from matrixInvw
%--------------------------------------------------------------------------------------------------------------------------
% written by Junzhong Xu, 2020-01-29
% 
%--------------------------------------------------------------------------------------------------------------------------

    %% preliminary
    vdist = cell([14 1]) ; 
    
    % ------------------ step#1: fitting cell-volumn-weighted cell size distribution -----------------
    % fitting
    [vdistAll, resnorm] = fit_lsqnonneg(signal, this.fitopts.Dictionary, this.fitopts.penalty) ; 
    % assign outputs
    matrixInvw = reshape(vdistAll(this.fitopts.indIn), [this.fitopts.NDin, this.fitopts.Nd]) ;    % intra
    matrixEx = reshape(vdistAll(this.fitopts.indEx),[this.fitopts.Nbetaex, this.fitopts.NDex]) ;      % extra
    vdist{1} = sum(matrixInvw,1) ; vdist{1} = vdist{1}' ; % d
    vdist{2} = sum(matrixInvw,2) ; % Din
    vdist{3} = sum(matrixEx,1) ; vdist{3} = vdist{3}' ; %Dex
    vdist{4} = sum(matrixEx,2) ; % betaex

    % ---------------------- step#2: fitting non-cell-volume-weighted cell size distribution --------------------
    % calculate intracellular signal
    Sin = this.fitopts.Dictionary(:,this.fitopts.indIn)*(vdistAll(this.fitopts.indIn)) ; vin = sum(matrixInvw(:)) ; 
    % fitting
    vdistAll2 = fit_lsqnonneg(Sin, this.fitopts.Dictionary2, this.fitopts.penalty2) ; 
    % apparent volumen fraction
    vdistAll2 = vdistAll2/sum(vdistAll2)*vin ; 
    % assign outputs
    matrixIn = reshape(vdistAll2(this.fitopts.indIn), [this.fitopts.NDin, this.fitopts.Nd]) ;    % intra
    vdist{5} = sum(matrixIn,1) ; vdist{5} = vdist{5}' ; % d
    vdist{6} = sum(matrixIn,2) ; % Din
    
    % ------------------------------  step#2, calculate dcal, Dincal divided by d^3 -------------------
    matrixInCal = matrixInvw ./ repmat(this.fitopts.ds'.^this.model.structure.Ndim,[this.fitopts.NDin,1]) ; 
    matrixInCal = matrixInCal / sum(matrixInCal(:)) * vin ; 
    vdist{7} = sum(matrixInCal,1) ; vdist{7} = vdist{7}' ; % dCal
    vdist{8} = sum(matrixInCal,2) ; % DinCal 
    
    % vdist{9} and vdist{10} TBD
    
    % ----------------------------------- output fitted matrices ----------------------------------------
    vdist{11} = matrixInvw ; vdist{12} = matrixEx ; vdist{13} = matrixIn ; vdist{14} = matrixInCal ; 

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vdist, resnorm] = fit_lsqnonneg(signal, Dictionary, Penalty)
% perform non-negative least-square fitting
%
% -------------------------------------------------------------------------------------------------------------------------

    % preparation
    H = Penalty*eye(size(Dictionary,2));
    d0 = zeros(size(Dictionary,2),1);
    temp=nonzeros(signal);
    % fitting
    [vdist, resnorm] = lsqnonneg([Dictionary;H],[temp;d0]) ; 

end
