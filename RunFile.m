addpath('Functions')
addpath('Functions/sMinimize')
%%
data='CM'
sz=1258;                % if you want to cross validate, pick less than the MSA size
%% Inference parameters
options.maxIter=500;

options.maxCycle=options.maxIter;     % for N=(MSA size) use 300, if using N=(MSA size/10) and lower use 800+
options.TolX=10^-5;     % how accurate do you need the parameters to be?
options.N=50;          % size of generated model MSA per BM iteration  
options.delta_t=100000; % burn time for each sequence
options.theta=0.3;      % reweighting for similarity, use 0.2-0.3
options.m=20;           % lbfgs steps to keep. the more the better, at the expense a bit of memory and time.
options.skip=ceil(options.maxIter/100);        % for logging purposes, write down the result every SKIP number of iteration into "output.wt"
                          % wt is an important log to check for convergence as a function of iter number
lambdaJ=0.01;             % coupling regularization
lambdah=0.01;             % field regularization
s_averaged=1;             % avg over stochstic runs. The lower the N, the more averaging you need.

%% Seed for CV set assignment
seed=4;
rng(seed)

%% Load align and divide to training and CV sets
load(['_Data/'   data '.mat'])

align=align(randperm(size(align,1)),:);
align_cv=align(sz+1:end,:);
align=align(1:sz,:);

%% Separate seed for inference
seed2=17
rng(seed2)


%% Initializations
q=max(max(align));
L=size(align,2);
J=zeros(q,q,L,L);
h=zeros(q,L);
N=size(align,1);

options.q=q;
options.L=L;

ws=[];

%% run the inference s_averaged times
for k=1:s_averaged
                        %   MSA,   regularization,  starting guess, 
    [J,h,N_eff,output]=SBM(align,  lambdaJ,lambdah,    J*0,h*0,    options); 
    ws(k,:)=Wj(J,h);
    
end

%% Finalize and save results
% average over the s_average runs
[J_av,h_av] = Jw(AvgOverRuns(ws,'discard')',q); % there is a throw away of outlying results. Maybe there shouldnt be?

% output structure is not being saved, 
% it takes a lot of space (because of wt). Only part is kept


% saved file will include a synthetic T=1 align
align_syn=Mex_IndepMonteCarlo(N,Wj(J_av,h_av),options.q,options.delta_t);

% and Trajectory summary for J and h separately
[J_tr,h_tr]=Trajectory(output.wt,q,0,options);
output_log=rmfield(output,'wt')

% save file containing all results and parameters (but no output w(t)!)
filename=[data '_m_' num2str(options.m) '_lam_' num2str(lambdaJ) '_Nmodel_' num2str(options.N) '_avg_' num2str(s_averaged) '_seed_' num2str(seed) '_L_' num2str(L) '_q_' num2str(q) '_N_' num2str(N) '_dt_' num2str(options.delta_t) '_iter_' num2str(options.maxIter) '.mat']
save(filename,'align','align_cv','align_syn','options','J_tr','h_tr','ws','J_av','h_av','s_averaged','lambdaJ','lambdah','output_log');

