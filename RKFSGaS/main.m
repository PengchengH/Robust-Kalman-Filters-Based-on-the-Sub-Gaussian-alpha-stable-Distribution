clc;clear;close all;
rng(10);

addpath('tools/dataAssimilation');
addpath('tools/tracking_model');
addpath('tools/stdt');
addpath('tools/SSG');
addpath('tools/gm');
addpath('tools/slash');
addpath('tools/vg');
addpath('tools/other');
addpath('tools/kftncm');
addpath('tools/dataset/modelling_result');
addpath('tools/dataset/traj');
addpath('tools/experiment_analysis');

% produce trajectory;
ratio_scale = 1;
Traj_produce;

% modelling index
% kftncm Kalman filter with true noise covariance matrices
modelling_index = {'kftncm','vg','slash','stdt','ssg'};
for i = 1:size(modelling_index,2)
    if ~exist(['result/' modelling_index{i}], 'dir')
        mkdir( ['result/' modelling_index{i}]);
    end
end

% Noise parameter set
noise = struct;
noise.model_index = {'gm','stdt','ssg'};
noise.dof_index   = [3   1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8;
                     0.3 0.5 0.7 0.9 1.2 1.7 2.5 3.5 6;
                     0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.85];


% Parameter index
particle_num_index = [10  40 100  400  1000  4000  10000];
point_num_index    = [1   2   4    6    10    20    30  ];
strategy_index     = {'GL','IS','IGGL','IGIS'};


num_base  = [size(noise.model_index,2),size(noise.dof_index,2),...
             size(modelling_index,2),size(particle_num_index,2),...
             size(point_num_index,2),size(strategy_index,2)];
num_index =  parallel_arrange(num_base);
           

% Parallel computation
n_proc = 4; 
c = parcluster('local');   % set clusters
c.NumWorkers = n_proc;
parpool(n_proc);
parfor ( i_index = 1:size(num_index,1), n_proc)
% for i_index = 1:size(num_index,1)
    Callfunction(num_index(i_index,:),modelling_index,noise,particle_num_index,...
        point_num_index,strategy_index);
end           

% delete parallel pool object
delete(gcp('nocreate'));

