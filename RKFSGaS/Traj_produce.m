rng(10);

addpath('tools/dataAssimilation');
addpath('tools/tracking_model');
addpath('tools/stdt');
addpath('tools/SSG');
addpath('tools/gm');
addpath('tools/slash');
addpath('tools/vg');
addpath('tools/other');
addpath('tools/dataset/modelling_result');
addpath('tools/dataset/traj');

% Experiemnt parameters
cycle = 100;

%% Produce model trajectory;
% modelling index
modelling_index = {'kftncm','vg','slash','stdt','ssg'};

% Noise parameter set
noise = struct;
noise.model_index = {'gm','stdt','ssg'};
noise.dof_index   = [3 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8;
                     0.3 0.5 0.7 0.9 1.2 1.7 2.5 3.5 6;
                     0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.85];

% Tracking model initialisation;
tracking_model_initial;

for noise_c =1:size(noise.model_index,2)
   for tail_c = 1 : size(noise.dof_index,2)  
       
       % Trajectory produce
       if strcmp(noise.model_index{noise_c},'gm')
          option    = struct;
          option.wp = 0.9;
          [traj] = tracking_model(noise.model_index{noise_c},...
                     noise.dof_index(1,tail_c),traj,cycle,option);
       elseif strcmp(noise.model_index{noise_c},'stdt')
          [traj] = tracking_model(noise.model_index{noise_c},...
                     noise.dof_index(2,tail_c),traj,cycle);
       elseif strcmp(noise.model_index{noise_c},'ssg')
          [traj] = tracking_model(noise.model_index{noise_c},...
                     noise.dof_index(3,tail_c),traj,cycle);
       end
       
       % Combine trajectory with modelling results
       for model_c = 1:size(modelling_index,2)
           
           if strcmp(modelling_index{model_c},'kftncm')
               save(['tools/dataset/traj/' 'md_' modelling_index{model_c}...
                      '_ns_' noise.model_index{noise_c}...
                      '_tail_' num2str(noise.dof_index(noise_c,tail_c)) '.mat'],'traj');               
           else
               load(['tools/dataset/modelling_result/' 'md_' modelling_index{model_c}...
                     '_ns_' noise.model_index{noise_c}...
                     '_tail_' num2str(noise.dof_index(noise_c,tail_c)) '.mat']); 

               md.Wscale = ratio_scale * md.Wscale;
               
               save(['tools/dataset/traj/' 'md_' modelling_index{model_c}...
                     '_ns_' noise.model_index{noise_c}...
                     '_tail_' num2str(noise.dof_index(noise_c,tail_c)) '.mat'],'traj','md','ns');
           end
       end
          
   end
end
