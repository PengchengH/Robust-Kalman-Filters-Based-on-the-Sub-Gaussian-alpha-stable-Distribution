function [fout] = post_process(fresult,traj)

if strcmp(fresult.state, 'success')
    
    fout  = struct;
    fout.state = 'success';

    x_error = fresult.x_pos-traj.x;

    prmse =  sqrt(mean(x_error(1,:,:).^2 + x_error(2,:,:).^2,'all'));
    vrmse =  sqrt(mean(x_error(3,:,:).^2 + x_error(4,:,:).^2,'all'));

    fout.prmse = prmse;
    fout.vrmse = vrmse;
    fout.time  = fresult.time;

    if isfield(fresult,'iteration') 
        ite_mean = mean(fresult.iteration(2:size(fresult.iteration,2),:),'all');
        fout.ite_mean = ite_mean;
        if isfield(fresult,'num_second')
            sed_ratio_mean = mean(fresult.num_second(2:size(fresult.iteration,2),:)./ ...
                    fresult.iteration(2:size(fresult.iteration,2),:),'all');
            fout.sed_ratio_mean = sed_ratio_mean;
        end
    end
    
elseif strcmp(fresult.state, 'fail')
    
    fout  = struct;
    fout.state = 'fail';
    
else
    
    error('the filtering result should include the right state');  
    
end


end

