function out = Fit(this, data)
%% mig.FitPars.Fit to fit any model to data
%
% -------------------------------------------------------------------------------------------------------------------------
% INPUTS:
%       this: a mig.FitPars object
%       data: a mig.ImageData object
% OUTPUTS:
%       out: a structure
%           .mapLabels: the name labels of output parametric map
%           .NumParms: number of fitting parameters
%           .{parm}: if fittingMethod=={'normal','general'}, fitted parameter maps. # of maps vary based on this.model.parms_labels
%           . resnorm: map of squared 2-norm of the residual 
%           .exitflag: map of exitflag. false: fitting successful; true: fitting failed
% -------------------------------------------------------------------------------------------------------------------------

    %% preliminary
    warning('off') ; 
    % check input 
    if ~isa(this,'mig.FitPars') || ~isa(data,'mig.ImageData'), error('%s: The inputs should be an FitPars object and an ImageData object',mfilename) ; end

    %% preprocessing data
    if this.fitopts.flag.denoise == 'y', end
    if this.fitopts.flag.degibbs == 'y', end
    if this.fitopts.flag.deivim == 'y', data = this.RemoveIVIM(data) ; disp('------ deivim performed -----------') ; end

    %% fitting preparation 
     switch this.fitopts.fittingMethod
        case {'normal','general'}       % normal/general fitting methods, i.e., minimize ObjectiveFcn     
            map0 = zeros([data.Nx, data.Ny, data.Nz]) ; 
            out.mapLabels = this.model.parmsLabels ; 
            for n = 1:length(out.mapLabels) , out.(out.mapLabels(n)) = map0 ; end
         case 'special'
            out.mapLabels = this.model.parmsLabels ; 
            out.parms = cell([data.Nx, data.Ny, data.Nz]) ; 
        otherwise
            error('%s: The fitopts.fittingMethod must be ''normal'' or ''special''',mfilename) ; 
     end
    out.resnorm = zeros([data.Nx, data.Ny, data.Nz]) ; out.exitflag = false([data.Nx, data.Ny, data.Nz]) ; 
    
    % options
    if this.fitopts.flag.parfor == 'y'
        poolobj = gcp ;  Ncore = poolobj.NumWorkers ; 
    else
        Ncore = 0 ; 
    end
    if this.fitopts.flag.multistart == 'y'  
        NumStarts = this.fitopts.NumStarts ; 
    else
        NumStarts = 1 ; 
    end
        
    %% fitting
    % calculate voxels with the mask
    [row, col, slice]=ind2sub(size(data.mask),find(data.mask>0)) ;     Npixel = length(row) ; 
    parmap1d = cell([Npixel,1]) ;     resnormmap1d = zeros([Npixel,1]) ; exitflag1d = zeros([Npixel,1]) ; 
    % prepare parfor
    if this.fitopts.flag.parfor == 'y'  % copy data and this(fitpars) to all workers to save communication time
        parData = parallel.pool.Constant(data) ; 
        parFitPars = parallel.pool.Constant(this) ; 
    else    % if no parfor, define variables for code consistency
        parData.Value = data ; 
        parFitPars.Value = this ; 
    end
    % ---------------------------------- fitting ------------------------------------ 
    parfor (npixel = 1:Npixel, Ncore)
%     for npixel=1:Npixel         % debugging purpose
        nx = row(npixel) ; ny = col(npixel) ; nz = slice(npixel) ; 
        signal = squeeze(parData.Value.img(nx,ny,nz,:)) ; 
        % noise sigma
        if isempty(parData.Value.sigma) || isscalar(parData.Value.sigma)
            sigma = parData.Value.sigma ; 
        else
            sigma = parData.Value.sigma(nx,ny,nz) ; 
        end
        % fitting
        try
           [x, resnorm] = run_optimization(parFitPars.Value, signal, sigma, NumStarts) ; 
        catch
            warning('Fitting in the pixel(%d,%d,%d) went wrong. Set all zeros for output parameters ...', nx, ny, nz)
            x = zeros([parFitPars.Value.model.NumParms,1]) ; resnorm = 0 ; 
            exitflag1d(npixel) = true ; 
        end
        % save fitting results
        parmap1d{npixel}  = x ; resnormmap1d(npixel) = resnorm ; 
    end % end of parfor loop
    % -------------------------------- output results ---------------------------------------------
    for npixel = 1:Npixel
        switch this.fitopts.fittingMethod
            case {'normal','general'}       % normal/general fitting methods, i.e., minimize ObjectiveFcn 
                for n=1:length(out.mapLabels)
                    out.(out.mapLabels(n))(row(npixel),col(npixel),slice(npixel)) = parmap1d{npixel}(n) ; 
                end
            case 'special'
                out.parms{row(npixel),col(npixel),slice(npixel)} = parmap1d{npixel} ; 
        end
        out.resnorm(row(npixel),col(npixel),slice(npixel)) = resnormmap1d(npixel) ; 
        out.exitflag(row(npixel),col(npixel),slice(npixel)) = exitflag1d(npixel) ; 
    end
    
end  % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [x, resnorm] = run_optimization(this,signal, sigma, NumStarts) 

    switch this.fitopts.fittingMethod
        case {'normal','general'}       % normal/general fitting methods, i.e., minimize ObjectiveFcn 
            % create optimization problem
            switch this.fitopts.solverName
                case {'lsqcurvefit'}         % the simple least-square curve fitting. It actualy calls lsqnonlin() 
                    problem = createOptimProblem('lsqcurvefit','x0',this.fitopts.x0,'objective',this.model.FcnSignal,'lb',this.fitopts.lb,'ub',this.fitopts.ub,'xdata',this.model,'ydata',signal,'options',this.fitopts.options);
                case {'lsqnonlin'}
                    ObjectiveFcn = @(parms)mig.FitPars.ObjectiveFcn(parms, this, signal, sigma) ; 
                    problem = createOptimProblem(this.fitopts.solverName,'x0',this.fitopts.x0,'objective',ObjectiveFcn,'lb',this.fitopts.lb,'ub',this.fitopts.ub,'options',this.fitopts.options);
                case {'fmincon'}
                    ObjectiveFcn = @(parms)mig.FitPars.ObjectiveFcn(parms, this, signal, sigma) ; 
                    problem = createOptimProblem(this.fitopts.solverName,'x0',this.fitopts.x0,'objective',ObjectiveFcn,'lb',this.fitopts.lb,'ub',this.fitopts.ub,'Aeq',this.fitopts.Aeq,'beq',this.fitopts.beq,'Aineq',this.fitopts.Aineq,'beq',this.fitopts.beq,'nonlcon',this.fitopts.nonlcon,'options',this.fitopts.options);
            end
            ms = MultiStart('UseParallel',false,'Display','off') ; 
            [x, resnorm] = run(ms,problem,NumStarts) ;    
        case 'special'          % special fitting method required by individual models
            [x, resnorm] = feval(this.fitopts.FitFcn, this, signal, sigma) ; 
    end
end

