function out = Fit(this, data)
%% mati.FitPars.Fit to fit any model to data
%
% -------------------------------------------------------------------------------------------------------------------------
% INPUTS:
%       this: a mati.FitPars object
%       data: a mati.ImageData object
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
    if ~isa(this,'mati.FitPars') || ~isa(data,'mati.ImageData'), error('%s: The inputs should be an FitPars object and an ImageData object',mfilename) ; end

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
     num_par=length(out.mapLabels);
    out.resnorm = zeros([data.Nx, data.Ny, data.Nz]) ; out.exitflag = false([data.Nx, data.Ny, data.Nz]) ; 
    out.MLI = zeros([data.Nx, data.Ny, data.Nz]) ;    % maximum likelihood
    out.minflag = zeros([data.Nx, data.Ny, data.Nz]) ; % whether the optimizaion is successful
    out.ncov = zeros([data.Nx, data.Ny, data.Nz, num_par, num_par]) ; % normalized covariance matrix 
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
    MLI1d = zeros([Npixel, 1]);    minflag1d = zeros([Npixel,1]) ;    ncov1d = cell([Npixel,1]) ;
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
    %  for npixel=1:Npixel         % debugging purpose
        nx = row(npixel) ; ny = col(npixel) ; nz = slice(npixel) ; 
        [nx ny nz]
        signal = squeeze(parData.Value.img(nx,ny,nz,:)) ; 
        % noise sigma
        if isempty(parData.Value.sigma) || isscalar(parData.Value.sigma)
            sigma = parData.Value.sigma ; 
        else
            sigma = parData.Value.sigma(nx,ny,nz) ; 
        end
        % fitting
        if isempty(nonzeros(signal))
            x = zeros([parFitPars.Value.model.NumParms,1]) ; resnorm = 0 ; MLI = 0; minflag = 0; ncov=0;
            exitflag1d(npixel) = true ;            
        else
            if (nx==2&&ny==54&&nz==12)
                sprintf('stop')
            end
            try
                [x, resnorm, MLI, minflag, ncov] = run_optimization(parFitPars.Value, signal, sigma, NumStarts) ; 
            catch error
                warning('Fitting in the pixel(%d,%d,%d) went wrong. Set all zeros for output parameters ...', nx, ny, nz)
                x = zeros([parFitPars.Value.model.NumParms,1]) ; resnorm = 0 ; MLI = 0; minflag = 0; ncov=0;
                exitflag1d(npixel) = true ; 
            end
        end
        % save fitting results
        parmap1d{npixel}  = x ; resnormmap1d(npixel) = resnorm ; MLI1d(npixel)=MLI; minflag1d(npixel)=minflag;
        ncov1d{npixel}  = reshape(ncov',1,[]) ;
    end % end of parfor loop
    % -------------------------------- output results ---------------------------------------------
    for npixel = 1:Npixel
        switch this.fitopts.fittingMethod
            case {'normal','general'}       % normal/general fitting methods, i.e., minimize ObjectiveFcn 
                for n=1:num_par                    
                    out.(out.mapLabels(n))(row(npixel),col(npixel),slice(npixel)) = parmap1d{npixel}(n) ;
                    out.MLI(row(npixel),col(npixel),slice(npixel))=MLI1d(npixel);
                    out.minflag(row(npixel),col(npixel),slice(npixel))=minflag1d(npixel);
                    try
                        out.ncov(row(npixel),col(npixel),slice(npixel),n,1:num_par) = ncov1d{npixel}((n-1)*num_par+1:n*num_par) ;                        
                    catch error
                        out.ncov(row(npixel),col(npixel),slice(npixel),n,1:num_par) = zeros([1,num_par]) ;
                        out.minflag(row(npixel),col(npixel),slice(npixel))=3;
                    end                   
                end
            case 'special'
                out.parms{row(npixel),col(npixel),slice(npixel)} = parmap1d{npixel} ; 
        end
        out.resnorm(row(npixel),col(npixel),slice(npixel)) = resnormmap1d(npixel) ; 
        out.exitflag(row(npixel),col(npixel),slice(npixel)) = exitflag1d(npixel) ; 
    end
    
end  % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [x, resnorm, MLI, minflag, ncov] = run_optimization(this,signal, sigma, NumStarts) 

    switch this.fitopts.fittingMethod
        case {'normal','general'}       % normal/general fitting methods, i.e., minimize ObjectiveFcn 
            % create optimization problem
            switch this.fitopts.solverName
                case {'lsqcurvefit'}         % the simple least-square curve fitting. It actualy calls lsqnonlin() 
                    problem = createOptimProblem('lsqcurvefit','x0',this.fitopts.x0,'objective',this.model.FcnSignal,'lb',this.fitopts.lb,'ub',this.fitopts.ub,'xdata',this.model,'ydata',signal,'options',this.fitopts.options);
                case {'lsqnonlin'}
                    ObjectiveFcn = @(parms)mati.FitPars.ObjectiveFcn(parms, this, signal, sigma) ; 
                    problem = createOptimProblem(this.fitopts.solverName,'x0',this.fitopts.x0,'objective',ObjectiveFcn,'lb',this.fitopts.lb,'ub',this.fitopts.ub,'options',this.fitopts.options);
                case {'fmincon'}
                    ObjectiveFcn = @(parms)mati.FitPars.ObjectiveFcn(parms, this, signal, sigma) ; 
                    problem = createOptimProblem(this.fitopts.solverName,'x0',this.fitopts.x0,'objective',ObjectiveFcn,'lb',this.fitopts.lb,'ub',this.fitopts.ub,'Aeq',this.fitopts.Aeq,'beq',this.fitopts.beq,'Aineq',this.fitopts.Aineq,'beq',this.fitopts.beq,'nonlcon',this.fitopts.nonlcon,'options',this.fitopts.options);
            end
            ms = MultiStart('UseParallel',false,'Display','off') ;
            [x,resnorm] = run(ms,problem,NumStarts) ;
            % we need the hessian matrix returned from fmincon to calculate
            % covariance matrix. However, the hessian matrix cannot be
            % returned when 'multi-start'option is used. So we repeat the
            % fmincon optimization one more time with a starting point
            % close to the fitted point from previous fmincon optimization.
            
            x(2)=x(2)+0.2; 
            [x,resnorm,~,~,~,~,hessian] = fmincon(ObjectiveFcn,x,this.fitopts.Aineq,this.fitopts.bineq,this.fitopts.Aeq,this.fitopts.beq,this.fitopts.lb,this.fitopts.ub,this.fitopts.nonlcon,this.fitopts.options);           
            [cov,ncov,minflag]=mati.Physics.hessianToCov(hessian,resnorm,length(signal),length(x));
            fval= this.model.FcnSignal(x, this.model);
            Lmax = mati.Physics.Get_Likelihood(signal,fval,sigma);
            MLI = Lmax*sqrt(abs(det(ncov)))*((2*pi).^(length(x)/2)); 
            MLI = MLI/prod(this.fitopts.ub-this.fitopts.lb);
        case 'special'          % special fitting method required by individual models
            [x, resnorm] = feval(this.fitopts.FitFcn, this, signal, sigma) ; 
    end
end

