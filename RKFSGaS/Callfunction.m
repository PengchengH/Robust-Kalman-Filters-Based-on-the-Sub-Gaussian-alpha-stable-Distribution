function [] = Callfunction(num_vry,modelling_index,noise,particle_num_index,...
        point_num_index,strategy_index)

    noise_c    = num_vry(1);
    tail_c     = num_vry(2);    
    model_c    = num_vry(3);
    particle_c = num_vry(4);
    point_c    = num_vry(5);
    strategy_c = num_vry(6);

    if strcmp(modelling_index{model_c},'ssg') ...
            &&  strcmp(noise.model_index{noise_c},'ssg')
  
        if strcmp(modelling_index{model_c},'ssg') 
            if (strcmp(strategy_index{strategy_c},'GL')  ...
                   || strcmp(strategy_index{strategy_c},'IGGL'))...   
                && particle_c == 1

                option = struct;
                option.optname=strategy_index{strategy_c};
                option.num_point= point_num_index(point_c);

                Callfunction_further(modelling_index{model_c},noise.model_index{noise_c},...
                noise.dof_index(noise_c,tail_c),option); 

            elseif (strcmp(strategy_index{strategy_c},'IS')  ...
                       || strcmp(strategy_index{strategy_c},'IGIS'))...   
                   && point_c == 1 

                option = struct;
                option.optname=strategy_index{strategy_c};
                option.num_particle= particle_num_index(particle_c);  

                Callfunction_further(modelling_index{model_c},noise.model_index{noise_c},...
                noise.dof_index(noise_c,tail_c),option); 

            end


        elseif ( strcmp(modelling_index{model_c},'kftncm') || ...
                 strcmp(modelling_index{model_c},'stdt') || ...
                 strcmp(modelling_index{model_c},'slash') || ...
                 strcmp(modelling_index{model_c},'vg')) && ...
               ( particle_c ==1 && point_c ==1 && strategy_c==1)

               Callfunction_further(modelling_index{model_c},...
               noise.model_index{noise_c},noise.dof_index(noise_c,tail_c));

        end

    end

end


function [] = Callfunction_further(model,noise,tail,varargin)
     
   if strcmp(model,'kftncm')
       % read data
       load(['tools/dataset/traj/'  'md_' model '_ns_' noise...
         '_tail_' num2str(tail) '.mat'],'traj');
       fresult = kf_call(traj);     
   else
       
       % read data
       load(['tools/dataset/traj/'  'md_' model '_ns_' noise...
         '_tail_' num2str(tail) '.mat'],'md','traj');
     
       % Filter Design
       if strcmp(model,'vg')
           GSM  = vg;
       elseif strcmp(model,'slash')
           GSM  = slash;
       elseif strcmp(model,'stdt')
           GSM  = stdt;
       elseif strcmp(model,'ssg')
           %% SSG selection extraction
           GSM  = ssg;
           if ~isempty(varargin)
               if isstruct(varargin{1})
                   opti = varargin{1};
                   if isfield(opti,'optname')   
                       GSM.option.optname = opti.optname;

                       % select the opt strategy
                       if strcmp(opti.optname,'GL')        
                           if isfield(opti,'num_point')
                               GSM.option.num_point = opti.num_point;        
                           end
                       elseif strcmp(opti.optname,'IS')
                           if isfield(opti,'num_particle')
                               GSM.option.num_particle = opti.num_particle;        
                           end
                       elseif strcmp(opti.optname,'IGGL')
                           if isfield(opti,'num_series')
                                 GSM.option.num_series = opti.num_series;        
                           end
                           if isfield(opti,'num_point')
                               GSM.option.num_point = opti.num_point;        
                           end                           
                       elseif strcmp(opti.optname,'IGIS')
                           if isfield(opti,'num_series')
                               GSM.option.num_series = opti.num_series;        
                           end   
                           if isfield(opti,'num_particle')
                               GSM.option.num_particle = opti.num_particle;        
                           end                          
                       else
                           error('the optname is not right');
                       end

                   else
                       error('the optname should be pointed or do not use option');
                   end
               else
                   error('the option should be a struct');
               end    
           end
           %%

       end
       [fresult] = RKF_GSM(md,traj,GSM);
       
   end
   
   [filterout] = post_process(fresult,traj);
   
   % save the results
   if strcmp(model,'ssg')
       
       if strcmp(GSM.option.optname,'GL') || strcmp(GSM.option.optname,'IGGL')
           save(['result/' model '/ns_' noise '_tail_' num2str(tail) '_'...
           GSM.option.optname '_point_' num2str(GSM.option.num_point) '.mat'],...
           'filterout');
       elseif strcmp(GSM.option.optname,'IS') || strcmp(GSM.option.optname,'IGIS')
           save(['result/' model '/noise_' noise '_tail_' num2str(tail) '_' ...
           GSM.option.optname '_npar_' num2str(GSM.option.num_particle) '.mat'],...
           'filterout');
       else
           error('GSM option error')
       end 
       
   else
       
        save(['result/' model '/ns_' noise '_tail_' num2str(tail) '.mat'],'filterout'); 
        
   end

end

