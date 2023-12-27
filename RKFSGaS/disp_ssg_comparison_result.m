clc;
clear;
close all;

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


% Noise parameter set
noise = struct;
noise.model_index = {'gm','stdt','ssg'};
noise.dof_index   = [3   1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8;
                     0.3 0.5 0.7 0.9 1.2 1.7 2.5 3.5 6;
                     0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7 1.85];

% Parameter index
particle_num_index = 100;
point_num_index    = 2;
strategy_index     = {'IS','GL','IGGL','IGIS'};


%% Analysis four SSG strategies
ssg_RMSE = cell(1,size(strategy_index,2));

for strategy_c = 1:size(strategy_index,2)
    
    % Extract data
    option = struct;
    option.optname = strategy_index{strategy_c};
    
    if strcmp(strategy_index{strategy_c},'GL') || strcmp(strategy_index{strategy_c},'IGGL')
        for point_c = 1:size(point_num_index,2)         
            option.num_point = point_num_index(point_c);
            ssg_RMSE{strategy_c} = RMSE_extract('ssg','ssg',noise.dof_index(3,:),...
                       'position','velocity','iteration','sec_ratio','time',option);
        end
        
    elseif strcmp(strategy_index{strategy_c},'IS') || strcmp(strategy_index{strategy_c},'IGIS')        
        for particle_c = 1:size(particle_num_index,2)         
            option.num_particle = particle_num_index(particle_c);
            ssg_RMSE{strategy_c} = RMSE_extract('ssg','ssg',noise.dof_index(3,:),...
                      'position','velocity','iteration','sec_ratio','time',option); 
        end  
    end
end

figure('Position',[50 50 500 500]);
tcl = tiledlayout(2,2);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';

nexttile;
hold on;
for strategy_c = 1:size(strategy_index,2)  
     plot(strategy_c*ones(2,1),ssg_RMSE{strategy_c}.scale,'--square',...
         'LineWidth',3,'MarkerSize',5);
end
hold off;
axis([1 size(strategy_index,2) 0 2]);
xlabel('(a)');
title('effective range','FontSize',14,'FontWeight','normal'); 

nexttile;
hold on;
for strategy_c = 1:size(strategy_index,2) 
     plot(ssg_RMSE{strategy_c}.prmse(1,:),ssg_RMSE{strategy_c}.prmse(2,:),'--square',...
         'LineWidth',3);
end
hold off;
xlabel('(b)');
title('position','FontSize',14,'FontWeight','normal');

nexttile;
hold on;
for strategy_c = 1:size(strategy_index,2) 
     plot(ssg_RMSE{strategy_c}.time(1,:),log(ssg_RMSE{strategy_c}.time(2,:)),'--square',...
         'LineWidth',3);
end
hold off;
xlabel('(c)');
title('implementation time','FontSize',14,'FontWeight','normal')

legend('IS','GL','IGGL','IGIS','FontSize',12,'Location',[0.65 0.2 0.22 0.2])

imaname=['RKFSGS_strategies' '.fig'];
saveas(gcf, imaname);  

imaname=['RKFSGS_strategies' '.png'];
saveas(gcf, imaname); 







