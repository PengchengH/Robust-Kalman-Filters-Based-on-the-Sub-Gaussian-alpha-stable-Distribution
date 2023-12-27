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
particle_num_index = [10  40 100  400  1000  4000  10000];
point_num_index    = [1   2   4    6    10    20    30  ];
strategy_index     = {'GL','IS','IGGL','IGIS'};


%% Analysis four SSG strategies
for strategy_c = 1:size(strategy_index,2)
    
    % Extract data
    option = struct;
    option.optname = strategy_index{strategy_c};
    
    if strcmp(strategy_index{strategy_c},'GL') || ...
        strcmp(strategy_index{strategy_c},'IGGL')
    
        ssg_RMSE = cell(1,size(point_num_index,2));
        
        for point_c = 1:size(point_num_index,2)         
            option.num_point = point_num_index(point_c);
            ssg_RMSE{point_c} = RMSE_extract('ssg','ssg',noise.dof_index(3,:),...
                 'position','velocity','iteration','sec_ratio','time',option);
        end
        
    elseif strcmp(strategy_index{strategy_c},'IS') || ...
        strcmp(strategy_index{strategy_c},'IGIS')
    
        ssg_RMSE = cell(1,size(particle_num_index,2));
        
        for particle_c = 1:size(particle_num_index,2)         
            option.num_particle = particle_num_index(particle_c);
            ssg_RMSE{particle_c} = RMSE_extract('ssg','ssg',noise.dof_index(3,:),...
                 'position','velocity','iteration','sec_ratio','time',option); 
        end  
    end
    
    % display data;
    if strcmp(strategy_index{strategy_c},'GL') || ...
        strcmp(strategy_index{strategy_c},'IGGL')
        if strcmp(strategy_index{strategy_c},'GL')
            figure('Position',[50 50 1000 750]);
        elseif strcmp(strategy_index{strategy_c},'IGGL')
            figure('Position',[50 50 1000 750+30]);
        end

        subplot(2,3,1);
        if strcmp(strategy_index{strategy_c},'GL')
            set (gca,'FontSize',14,'position',[0+0.06,1/2+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGGL')
            set (gca,'FontSize',14,'position',[0+0.06,0.03+1/2+0.15,0.25,0.285] );
        end

        hold on;
        for point_c = 1:size(point_num_index,2)  
             plot(point_c*ones(2,1),ssg_RMSE{point_c}.scale,'--square',...
                 'DisplayName',['L=' num2str(point_num_index(point_c))],...
                 'LineWidth',3,'MarkerSize',5);
        end
        hold off;
        axis([0 size(point_num_index,2) 0 2]);
        xlabel('point number index');
        ylabel('shape parameter');
        title('effective range'); 
        text(0.5,-0.35,'(a)','FontSize',14, 'Units','normalized')

        subplot(2,3,2);
        if strcmp(strategy_index{strategy_c},'GL')
            set (gca,'FontSize',14,'position',[1/3+0.06,1/2+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGGL')
            set (gca,'FontSize',14,'position',[1/3+0.06,0.03+1/2+0.15,0.25,0.285] );
        end

        hold on;
        for point_c = 1:size(point_num_index,2)  
             plot(ssg_RMSE{point_c}.prmse(1,:),ssg_RMSE{point_c}.prmse(2,:),'--',...
                 'DisplayName',['L=' num2str(point_num_index(point_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('RMSE');
        title('position');
        text(0.5,-0.35,'(b)','FontSize',14, 'Units','normalized')
        
        subplot(2,3,3);
        if strcmp(strategy_index{strategy_c},'GL')
            set (gca,'FontSize',14,'position',[2/3+0.06,1/2+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGGL')
            set (gca,'FontSize',14,'position',[2/3+0.06,0.03+1/2+0.15,0.25,0.285] );
        end
        hold on;
        for point_c = 1:size(point_num_index,2)  
             plot(ssg_RMSE{point_c}.vrmse(1,:),ssg_RMSE{point_c}.vrmse(2,:),'--',...
                 'DisplayName',['L=' num2str(point_num_index(point_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('RMSE');
        title('velocity');  
        text(0.5,-0.35,'(c)','FontSize',14, 'Units','normalized')
        
        subplot(2,3,4);
        if strcmp(strategy_index{strategy_c},'GL')
            set (gca,'FontSize',14,'position',[0+0.06,0+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGGL')
            set (gca,'FontSize',14,'position',[0+0.06,0.075+0+0.14,0.25,0.285] );
        end
        hold on;
        for point_c = 1:size(point_num_index,2)  
             plot(ssg_RMSE{point_c}.ite_mean(1,:),ssg_RMSE{point_c}.ite_mean(2,:),'--',...
                 'DisplayName',['L=' num2str(point_num_index(point_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('iteration number');
        title('average iteration number');
        text(0.5,-0.35,'(d)','FontSize',14, 'Units','normalized')

        
        subplot(2,3,5);
        if strcmp(strategy_index{strategy_c},'GL')
            set (gca,'FontSize',14,'position',[1/3+0.06,0+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGGL')
            set (gca,'FontSize',14,'position',[1/3+0.06,0.075+0+0.14,0.25,0.285] );
        end
        hold on;
        for point_c = 1:size(point_num_index,2)  
             plot(ssg_RMSE{point_c}.time(1,:),log(ssg_RMSE{point_c}.time(2,:)),'--',...
                 'DisplayName',['L=' num2str(point_num_index(point_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('log(Time)');
        title('implementation time'); 
        text(0.5,-0.35,'(e)','FontSize',14, 'Units','normalized')


        
        if strcmp(strategy_index{strategy_c},'IGGL')
            subplot(2,3,6)
            if strcmp(strategy_index{strategy_c},'GL')
                set (gca,'FontSize',14,'position',[2/3+0.06,0+0.15,0.25,0.315] );
            elseif strcmp(strategy_index{strategy_c},'IGGL')
                set (gca,'FontSize',14,'position',[2/3+0.06,0.075+0+0.14,0.25,0.285] );
            end
            hold on;
            for point_c = 1:size(point_num_index,2)  
                 plot(ssg_RMSE{point_c}.sed_ratio_mean(1,:),ssg_RMSE{point_c}.sed_ratio_mean(2,:),'--',...
                     'DisplayName',['L=' num2str(point_num_index(point_c))],'LineWidth',3);
            end
            hold off;
            xlabel('shape parameter');
            ylabel('ratio');
            title('second method ratio'); 
            text(0.5,-0.35,'(f)','FontSize',14, 'Units','normalized')
            
            lgd = legend('Location',[0.12 0.02 0.75 0.05]);
            lgd.FontSize = 12;
            lgd.FontWeight = 'bold';
            lgd.Orientation = 'horizontal';
        end
        
        if strcmp(strategy_index{strategy_c},'GL')
        
            lgd = legend('Location',[0.8 0.1 0.05 0.3]);
            lgd.FontSize = 12;
            lgd.FontWeight = 'bold';
            
        end
        
    end
    
%%
    if strcmp(strategy_index{strategy_c},'IS') || ...
        strcmp(strategy_index{strategy_c},'IGIS')

        if strcmp(strategy_index{strategy_c},'IS')
            figure('Position',[50 50 1000 750]);
        elseif strcmp(strategy_index{strategy_c},'IGIS')
            figure('Position',[50 50 1000 750+30]);
        end

        subplot(2,3,1);
        if strcmp(strategy_index{strategy_c},'IS')
            set (gca,'FontSize',14,'position',[0+0.06,1/2+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGIS')
            set (gca,'FontSize',14,'position',[0+0.06,0.03+1/2+0.15,0.25,0.285] );
        end 

        hold on;
        for particle_c = 1:size(particle_num_index,2)  
             plot(particle_c*ones(2,1),ssg_RMSE{particle_c}.scale,'--square',...
                 'DisplayName',['N=' num2str(particle_num_index(particle_c))],...
                 'LineWidth',3,'MarkerSize',5);
        end
        hold off;
        axis([0 size(particle_num_index,2) 0 2]);
        xlabel('particle number index');
        ylabel('shape parameter');
        title('effective range'); 
        text(0.5,-0.35,'(a)','FontSize',14, 'Units','normalized')

        subplot(2,3,2);
        if strcmp(strategy_index{strategy_c},'IS')
            set (gca,'FontSize',14,'position',[1/3+0.06,1/2+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGIS')
            set (gca,'FontSize',14,'position',[1/3+0.06,0.03+1/2+0.15,0.25,0.285] );
        end       
        hold on;
        for particle_c = 1:size(particle_num_index,2)  
             plot(ssg_RMSE{particle_c}.prmse(1,:),ssg_RMSE{particle_c}.prmse(2,:),'--',...
                 'DisplayName',['N=' num2str(particle_num_index(particle_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('RMSE');
        title('position');
        text(0.5,-0.35,'(b)','FontSize',14, 'Units','normalized')

        subplot(2,3,3);
        if strcmp(strategy_index{strategy_c},'IS')
            set (gca,'FontSize',14,'position',[2/3+0.06,1/2+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGIS')
            set (gca,'FontSize',14,'position',[2/3+0.06,0.03+1/2+0.15,0.25,0.285] );
        end

        hold on;
        for particle_c = 1:size(particle_num_index,2)  
             plot(ssg_RMSE{particle_c}.vrmse(1,:),ssg_RMSE{particle_c}.vrmse(2,:),'--',...
                 'DisplayName',['N=' num2str(particle_num_index(particle_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('RMSE');
        title('velocity');  
        text(0.5,-0.35,'(c)','FontSize',14, 'Units','normalized')
        
        subplot(2,3,4);
        if strcmp(strategy_index{strategy_c},'IS')
            set (gca,'FontSize',14,'position',[0+0.06,0+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGIS')
            set (gca,'FontSize',14,'position',[0+0.06,0.075+0+0.14,0.25,0.285] );
        end
        hold on;
        for particle_c = 1:size(particle_num_index,2)  
             plot(ssg_RMSE{particle_c}.ite_mean(1,:),ssg_RMSE{particle_c}.ite_mean(2,:),'--',...
                 'DisplayName',['N=' num2str(particle_num_index(particle_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('iteration number');
        title('average iteration number');
        text(0.5,-0.35,'(d)','FontSize',14, 'Units','normalized')
        
        
        subplot(2,3,5);
        if strcmp(strategy_index{strategy_c},'IS')
            set (gca,'FontSize',14,'position',[1/3+0.06,0+0.15,0.25,0.315] );
        elseif strcmp(strategy_index{strategy_c},'IGIS')
            set (gca,'FontSize',14,'position',[1/3+0.06,0.075+0+0.14,0.25,0.285] );
        end
        hold on;
        for particle_c = 1:size(particle_num_index,2)  
             plot(ssg_RMSE{particle_c}.time(1,:),log(ssg_RMSE{particle_c}.time(2,:)),'--',...
                 'DisplayName',['N=' num2str(particle_num_index(particle_c))],'LineWidth',3);
        end
        hold off;
        xlabel('shape parameter');
        ylabel('log(Time)');
        title('implementation time'); 
        text(0.5,-0.35,'(e)','FontSize',14, 'Units','normalized')

        if strcmp(strategy_index{strategy_c},'IGIS')
            subplot(2,3,6);
            if strcmp(strategy_index{strategy_c},'IS')
                set (gca,'FontSize',14,'position',[2/3+0.06,0+0.15,0.25,0.315] );
            elseif strcmp(strategy_index{strategy_c},'IGIS')
                set (gca,'FontSize',14,'position',[2/3+0.06,0.075+0+0.14,0.25,0.285] );
            end
            hold on;
            for particle_c = 1:size(particle_num_index,2)  
                 plot(ssg_RMSE{particle_c}.sed_ratio_mean(1,:),ssg_RMSE{particle_c}.sed_ratio_mean(2,:),'--',...
                     'DisplayName',['N=' num2str(particle_num_index(particle_c))],'LineWidth',3);
            end
            hold off;
            xlabel('shape parameter');
            ylabel('ratio');
            title('second method ratio'); 
            text(0.5,-0.35,'(f)','FontSize',14, 'Units','normalized')
            
            lgd = legend('Location',[0.12 0.02 0.75 0.05]);
            lgd.FontSize = 12;
            lgd.FontWeight = 'bold';
            lgd.Orientation = 'horizontal';
            lgd.NumColumns  = 4;
            
        end
        
        if strcmp(strategy_index{strategy_c},'IS')
        
            lgd = legend('Location',[0.8 0.1 0.05 0.3]);
            lgd.FontSize = 12;
            lgd.FontWeight = 'bold';
            
        end
        

        
    end    
 
 %%
    imaname=['RKFSGS_' strategy_index{strategy_c} '.fig'];
    saveas(gcf, imaname);  

    imaname=['RKFSGS_' strategy_index{strategy_c} '.png'];
    saveas(gcf, imaname);  
 
end






