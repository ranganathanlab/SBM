function [J,h,N_eff,output]=SBM(align,lambdaJ,lambdah,J_init,h_init,options)
    
    %similarity weighting
    if options.theta>0
        [W,N_eff]=Mex_CalcWeights(align,options.theta);
    else
        N_eff=size(align,1);
        W=zeros(size(align,1),1)+1;
    end
    fprintf(['N_{eff}=' num2str(N_eff) '\n'])
    
    %assign default size of model sample and default features
    if ~isfield(options,'N')
        options.N=N_eff;
    end
    if ~isfield(options,'prune')
        options.prune=1;
    end  
    if ~isfield(options,'regmat')
        options.regmat=0;
    end 
    if ~isfield(options,'sim_weight')
        options.sim_weight=0;
    end 
    

    %find data statistics
    [fi,fij]=CalcStatsWeighted(options.q,align,W/N_eff);   
    
    %Minimization procedure
    w0=Wj(J_init,h_init);
    [w,output]=sMinimizer(@(x) LogLike(x,lambdaJ,lambdah,fi,fij,options),w0,options);
    [J,h]=Jw(w,options.q);
    
    %gauge final result, just in case
    [J,h]=IsingGauge(J,h);
end


function [grad]=LogLike(wr,lambdaJ,lambdah,fi,fij,options)
    [J,h]=Jw(wr,options.q);

    align_mod=Mex_IndepMonteCarlo(options.N,Wj(J,h),options.q,options.delta_t);
    
    % now gradient of it
    if options.sim_weight==1
        [p,N_eff]=Mex_CalcWeights(align_mod,options.theta);
    else        
        p=zeros(options.N,1)+1/options.N;
    end
    %model MSA stats
    [fi_tot,fij_tot]=CalcStatsWeighted(options.q,align_mod,p);
    
    %compared to data MSA stats - that is the update rule
    gradh=fi_tot-fi+2*lambdah*h;
    gradJ=(fij_tot-fij+2*lambdaJ*J+2*options.regmat*J).*options.prune;
    grad=Wj(gradJ,gradh);
   
end
