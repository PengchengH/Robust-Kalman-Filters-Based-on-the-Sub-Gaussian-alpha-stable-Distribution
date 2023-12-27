function [] = callmodel(GMn,stdtn,SSGn,tail_c,model_name,noise_name,n,Cycle)

% select noise model
if strcmp(noise_name,'gm')
    Q = GMn.Q;
    d = size(Q,2);
    U = GMn.U_index(tail_c);
    p = GMn.p;
    
    Q1   = Q;
    Q2   = U*Q;
    
    noise      = struct;
    noise.Q1   = Q1;
    noise.Q2   = Q2;
    noise.U    = U;
    noise.name = noise_name;
    noise.samplesize = n;
    
elseif strcmp(noise_name,'stdt')
    Q   = stdtn.Q;
    d   = size(Q,2);
    dof = stdtn.dof_index(tail_c);
    
    noise      = struct;
    noise.Q    = Q;
    noise.dof  = dof;
    noise.name = noise_name;
    noise.samplesize = n;
elseif strcmp(noise_name,'ssg')
    Q   = SSGn.Q;
    d   = size(Q,2);
    dof = SSGn.dof_index(tail_c);
    
    noise      = struct;
    noise.Q    = Q;
    noise.dof  = dof;
    noise.name = noise_name;
    noise.samplesize = n;
end


% modelling fitting
if strcmp(model_name,'vg_EM')
    fitter = vg_modelling;
elseif strcmp(model_name,'vg_MLE')
    fitter = vg_modelling;
elseif strcmp(model_name,'vg_moment')
    fitter = vg_modelling;  
elseif strcmp(model_name,'slash_EM')
    fitter = slash_modelling;
elseif strcmp(model_name,'slash_MLE')
    fitter = slash_modelling;
elseif strcmp(model_name,'stdt_EM')
    fitter = stdt_modelling;
elseif strcmp(model_name,'stdt_MLE')
    fitter = stdt_modelling;
elseif strcmp(model_name,'SSG_ECF')
    fitter = SSG_modelling;
elseif strcmp(model_name,'SSG_EM')
    fitter = SSG_modelling;
elseif strcmp(model_name,'SSG_MLE')
    fitter = SSG_modelling;
elseif strcmp(model_name,'SSG_Mellin')
    fitter = SSG_modelling;
elseif strcmp(model_name,'SSG_Percentile')
    fitter = SSG_modelling;
end

% ESTIMATION PARAMATER
model       = struct;
model.name  = model_name;
model.SIGMA = zeros(d,d,Cycle);
model.dof   = zeros(Cycle,1);
model.time  = zeros(1,1);


tic; % timer
for Cycle_c = 1:Cycle
   % sampling
   if strcmp(noise_name,'gm')
       Samp = GM_noise(Q1,Q2,p,n);
   elseif strcmp(noise_name,'stdt')
       Samp = strnd(zeros(1,d),Q,dof,n)';
   elseif strcmp(noise_name,'ssg')
       Samp = SSGrdn(zeros(d,1),Q,dof,n);
   end
   
   A = chol(Q)';
   Y = reshape(inv(A)*Samp,1,[]);
   
   para_1 = fitter.fitting(Y,model_name);
   model.SIGMA(:,:,Cycle_c) = para_1.SIGMA*Q;
   model.dof(Cycle_c)       = para_1.dof;
%    disp(Cycle_c);
end
model.time = toc;

%   figure(2);
%   subplot(1,2,1);
%   plot(model.dof);
%   subplot(1,2,2);
%   plot(squeeze(model.SIGMA(1,1,:)));


save(['result/' model_name '_' noise_name '_tail_' num2str(tail_c) '.mat'],'model','noise');
end

