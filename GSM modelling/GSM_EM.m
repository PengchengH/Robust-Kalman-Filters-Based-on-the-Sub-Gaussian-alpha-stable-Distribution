clc;
clear;
close all;
warning('off');

addpath('tools/other');
addpath('tools/SLASH');
addpath('tools/STDT');
addpath('tools/SSG');
addpath('tools/VG');
addpath('tools/GM');

%% Initialisation
T  = 1; 
Qn = 10*diag(ones(2,1)); 
n  = 1000; 

% noise parameter set
GMn    = struct;
GMn.p  = 0.9;
GMn.Q  = Qn;
GMn.U_index = [3 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8];

stdtn   = struct;
stdtn.Q = Qn;
stdtn.dof_index = [0.3 0.5 0.7 0.9 1.2 1.7 2.5 3.5 6];


SSGn    = struct;
SSGn.Q  = Qn;
SSGn.dof_index = [1.85];

num_tail = size(GMn.U_index,2);

Cycle = 1000;

noise_index = {'gm','stdt','ssg'};
model_index = {'vg_MLE','vg_moment','slash_MLE','stdt_MLE','SSG_EM','SSG_Mellin'};
% model_index = {'vg_MLE','vg_moment','slash_MLE','stdt_MLE','SSG_Mellin'};

% variable number set
num_vry  = [1,1,0];
num_base = [num_tail,size(model_index,2),size(noise_index,2)];
i_sum    = prod(num_base);

num_vry_index = zeros(i_sum,size(num_vry,2));
for i_index = 1:i_sum
    num_vry(size(num_vry,2)) = num_vry(size(num_vry,2))+1;
    for i = 1:size(num_vry,2)
        j = size(num_vry,2)+1-i;
        if num_vry(j) > num_base(j)
           num_vry(j)   = 1;
           num_vry(j-1) = num_vry(j-1)+1;
        end
    end
    num_vry_index(i_index,:) = num_vry;
end

n_proc =8;
c = parcluster('local');   % set clusters
c.NumWorkers = n_proc;
parpool(n_proc);
parfor (i_index = 1:i_sum, n_proc)
% for i_index = 1:i_sum
    num_vry = num_vry_index(i_index,:);
    tail_c  = num_vry(1);
    model_c = num_vry(2);
    noise_c   = num_vry(3);   
    
%       if model_c == 1  && noise_c == 3 && tail_c ==1
             callmodel(GMn,stdtn,SSGn, tail_c,model_index{model_c},...
                noise_index{noise_c},n,Cycle);
%       end
    
end

% delete parallel pool object
delete(gcp('nocreate'));













