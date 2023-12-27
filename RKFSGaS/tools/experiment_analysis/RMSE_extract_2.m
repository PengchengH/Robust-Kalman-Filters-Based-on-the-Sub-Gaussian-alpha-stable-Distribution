function [out] = RMSE_extract(model,noise,tail_index,varargin)

out = struct;

if strcmp(model,'ssg')
    if isstruct(varargin{end})
        opti = varargin{end};
        if isfield(opti,'optname') 
            if strcmp(opti.optname,'GL') || strcmp(opti.optname,'IGGL')        
                if ~isfield(opti,'num_point')
                    error('num_point should be inputed for GL and IGGL')
                end
            elseif strcmp(opti.optname,'IS') || strcmp(opti.optname,'IGIS') 
                if ~isfield(opti,'num_particle')
                    error('num_particle should be inputed for IS and IGIS')
                end 
            else
                error('the optname is illegal')
            end
            
        else
            error('the optname should be inputed');
        end
        
    else
        error('the option for ssg is necessary and should be struct');
    end
end

if strcmp(measure,'position')
    out.prmse = zeros(2,size(tail_index,2));
    out.prmse(1,:) = tail_index;
elseif  strcmp(measure,'velocity')
    out.vrmse = zeros(2,size(tail_index,2));
    out.vrmse(1,:) = tail_index;
elseif strcmp(measure,'iteration')
    out.ite_mean = zeros(2,size(tail_index,2));
    out.ite_mean(1,:) = tail_index;
elseif strcmp(measure,'sec_ratio')
    out.sed_ratio_mean = zeros(2,size(tail_index,2));
    out.sed_ratio_mean(1,:) = tail_index;        
end



for tail_c = 1:size(tail_index,2)
    
    tail = tail_index(tail_c);
    
    if strcmp(model,'ssg')
        if strcmp(opti.optname,'GL') || strcmp(opti.optname,'IGGL')
            load(['result/' model '/ns_' noise '_tail_' num2str(tail) '_'...
                 opti.optname '_point_' num2str(opti.num_point) '.mat']);
        elseif strcmp(opti.optname,'IS') || strcmp(opti.optname,'IGIS')
            load(['result/' model '/noise_' noise '_tail_' num2str(tail) '_' ...
                 opti.optname '_npar_' num2str(opti.num_particle) '.mat']);
        end    
    else
        load(['result/' model '/ns_' noise '_tail_' num2str(tail) '.mat']); 
    end
    
    if isfield(out,'prmse')  
        out.prmse(2,tail_c) = filterout.prmse;
    end
    if isfield(out,'vrmse')  
        out.vrmse(2,tail_c) = filterout.vrmse;
    end
    if isfield(out,'ite_mean')  
        out.ite_mean(2,tail_c) = filterout.ite_mean;
    end  
    if isfield(out,'sed_ratio_mean')  
        out.sed_ratio_mean(2,tail_c) = filterout.sed_ratio_mean;
    end      
end


end

