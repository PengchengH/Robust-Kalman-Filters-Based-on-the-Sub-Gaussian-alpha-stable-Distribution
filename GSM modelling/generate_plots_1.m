clc;
clear;
close all;
rng(10);

addpath('tools/other');
addpath('tools/SLASH');
addpath('tools/STDT');
addpath('tools/SSG');
addpath('tools/VG');
addpath('tools/GM');


%% Sampling initialisation
%% Initialisation
T  = 1; 
d  = 2;
Qn = 10*diag(ones(d,1)); 
n  = 1000; 

% noise parameter set
GMn    = struct;
GMn.p  = 0.9;
GMn.Q  = Qn;
GMn.U_index = [1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8];

stdtn   = struct;
stdtn.Q = Qn;
stdtn.dof_index = [0.3 0.5 0.7 0.9 1.2 1.7 2.5 3.5];

SSGn    = struct;
SSGn.Q  = Qn;
SSGn.dof_index = [0.3 0.5 0.7 0.9 1.1 1.3 1.5 1.7];

num_tail = 9;

Cycle = 1000;

noise_index = {'gm','stdt','ssg'};
model_index = {'vg_MLE','slash_MLE','stdt_MLE','SSG_EM'};
model_index_1 = {'vg','slash','stdt','SSG'};

dof_data_error  = zeros(num_tail,size(noise_index,2),size(model_index,2));
delta_mean        = zeros(num_tail,size(noise_index,2),size(model_index,2));
sigma_data_error  = zeros(d,d,num_tail,size(noise_index,2),size(model_index,2)); 
W_scale  = zeros(d,d,num_tail,size(noise_index,2),size(model_index,2));
W_dof  = zeros(num_tail,size(noise_index,2),size(model_index,2));


for tail_c = 1:num_tail
for noise_c = 1:size(noise_index,2)
for model_c = 1:size(model_index,2)
 

if isfile(['result/' model_index{model_c} '_' noise_index{noise_c} '_tail_' num2str(tail_c) ...
    '.mat'])

load(['result/' model_index{model_c} '_' noise_index{noise_c} '_tail_' num2str(tail_c) ...
       '.mat']);

    dof_data_error(tail_c,noise_c ,model_c)   = rms((model.dof - mean(model.dof))); 
    delta_mean(tail_c,noise_c ,model_c)       = mean(model.dof);
    sigma_data_error(tail_c,noise_c ,model_c) = rms((model.SIGMA - mean(model.SIGMA,3)),"all");
    
    SIGMA = model.SIGMA;
    for i = 1:size(model.SIGMA,3)
        SIGMA(:,:,i) = inv(model.SIGMA(:,:,i));
    end
    
    [W_dof(tail_c,noise_c ,model_c),W_scale(:,:,tail_c,noise_c ,model_c)]...
           = Wishart_modelling(SIGMA);
    W_scale(:,:,tail_c,noise_c ,model_c) = inv(W_scale(:,:,tail_c,noise_c ,model_c));
    
    md = struct;
    md.name   =model_index_1{model_c};
    md.modelling   = model_index{model_c};
    md.delta  = delta_mean(tail_c,noise_c ,model_c);
    md.Wdof   = W_dof(tail_c,noise_c ,model_c);
    md.Wscale = W_scale(:,:,tail_c,noise_c ,model_c);
    
    ns     = noise;
    
    if isfield(ns,'U')
        save(['result_17/' 'md_' model_index_1{model_c} '_ns_' ...
       noise_index{noise_c} '_tail_' num2str(ns.U) '.mat'],'md','ns');
    elseif isfield(ns,'dof')
        save(['result_17/' 'md_' model_index_1{model_c} '_ns_' ...
       noise_index{noise_c} '_tail_' num2str(ns.dof) '.mat'],'md','ns');
    else
        error('there is no dof value in the struct');
    end
       
else
    dof_data_error(tail_c,noise_c ,model_c)   = NaN;  
    sigma_data_error(tail_c,noise_c ,model_c) = NaN;
end


end
end
end

noise_index = {'GM','ST',['SG' '\alpha' 'S']};
model_index = {'VG','SL','ST',['SG' '\alpha' 'S']};


%% display
figure('Position',[100 50 1000 2/3*1000]);
for noise_c = 1:size(noise_index,2)
%     if strcmp(noise_index{noise_c},'GM')
%         x_index = GMn.U_index;
%     elseif strcmp(noise_index{noise_c},'stdt')
%         x_index = stdtn.dof_index;
%     elseif strcmp(noise_index{noise_c},'SSG')
%         x_index = SSGn.dof_index;
%     end
for model_c = 1:size(model_index,2)
    subplot(size(noise_index,2),size(model_index,2),(noise_c-1)*size(model_index,2)+model_c);
    semilogy(1:1:num_tail,abs(dof_data_error(:,noise_c,model_c)),'-','LineWidth',3);
    hold on;
    semilogy(1:1:num_tail,abs(delta_mean(:,noise_c,model_c)),'-','LineWidth',3);
    hold off;
    title([noise_index{noise_c} '&' model_index{model_c}]);
end  
end
legend('RMSE','MEAN','Orientation','horizontal','FontSize',12,'Location',[0.4 0.03 0.2 0]);
imaname=['delta-dof-modelling' '.png'];
saveas(gcf, imaname);  
imaname=['delta-dof-modelling' '.fig'];
saveas(gcf, imaname); 
% figure(2);
% for noise_c = 1:size(noise_index,2)
% for model_c = 1:size(model_index,2)
%     subplot(size(noise_index,2),size(model_index,2),(noise_c-1)*size(model_index,2)+model_c);
%     semilogy(1:1:num_tail,abs(sigma_data_error(:,noise_c,model_c)),'-ob','LineWidth',2);
%     title(['Sigma-rms-' noise_index{noise_c} '&' model_index{model_c}]);
% end  
% end

% figure(3);
% for noise_c = 1:size(noise_index,2)
% for model_c = 1:size(model_index,2)
%     subplot(size(noise_index,2),size(model_index,2),(noise_c-1)*size(model_index,2)+model_c);
%     semilogy(1:1:num_tail,abs(delta_mean(:,noise_c,model_c)),'-or','LineWidth',2);
%     title(['delta-mean-' noise_index{noise_c} '&' model_index{model_c}]);
% end  
% end

figure('Position',[100 50 1000 2/3*1000]);
for noise_c = 1:size(noise_index,2)
for model_c = 1:size(model_index,2)
    subplot(size(noise_index,2),size(model_index,2),(noise_c-1)*size(model_index,2)+model_c);
    semilogy(1:1:num_tail,abs(W_dof(:,noise_c,model_c)),'-',"Color","#D95319",'LineWidth',3);
    title([noise_index{noise_c} '&' model_index{model_c}]);
end  
end
imaname=['IW-dof-modelling' '.png'];
saveas(gcf, imaname); 
imaname=['IW-dof-modelling' '.fig'];
saveas(gcf, imaname); 

figure('Position',[100 50 1000 2/3*1000]);
for noise_c = 1:size(noise_index,2)
for model_c = 1:size(model_index,2)
    subplot(size(noise_index,2),size(model_index,2),(noise_c-1)*size(model_index,2)+model_c);
    semilogy(1:1:num_tail,abs(squeeze(W_scale(1,1,:,noise_c,model_c))),'-',"Color","#D95319",'LineWidth',3);
    title([noise_index{noise_c} '&' model_index{model_c}]);
end  
end
imaname=['IW-sclae-modelling' '.png'];
saveas(gcf, imaname); 
imaname=['IW-sclae-modelling' '.fig'];
saveas(gcf, imaname); 



% -- END OF FILE --
