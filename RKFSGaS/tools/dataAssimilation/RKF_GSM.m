function [RKF] = RKF_GSM(para,traj,GSM)

% Parameter reading
n = size(traj.Qs,1);
% m = size(para.SIGMA,1);
length = traj.length;
cycle = traj.cycle;

% GSM_filter related matrix
x_pos = zeros(n,length,cycle);
P_pos = zeros(n,n,length,cycle);

for j = 1:cycle
    x_pos(:,1,j)   = traj.x0(:,j);
    P_pos(:,:,1,j) = traj.P0(:,:,j);
end

K_lamda = zeros(length,cycle);
ite     = zeros(length,cycle);
num_second = zeros(length,cycle);

state  = 'success';

tic; % timer
for cycle_c = 1:cycle 
    for t = 2:traj.length
        % extract parameter
        mea =  traj.z(:,t,cycle_c);
        xt_ =  x_pos(:,t-1,cycle_c);
        Pt_ =  P_pos(:,:,t-1,cycle_c);

        [out] = assimation(mea,xt_,Pt_,para,traj,GSM);

        x_pos(:,t,cycle_c)  = out.x_pos;
        P_pos(:,:,t,cycle_c) = out.P_pos;
        K_lamda(t,cycle_c)   = out.K_lamda;
        ite(t,cycle_c)       = out.ite;
        num_second(t,cycle_c) = out.num_second;

%         disp(out.P_pos);
        
        if strcmp(out.state, 'fail')
            state = out.state;
            break;
        end
%          disp(t);
    end
    
    if strcmp(state, 'fail')
        break;
    end

% figure(1);
% subplot(2,2,1);
% plot(traj.x(1,:,cycle_c),traj.x(2,:,cycle_c),'-*k');
% hold on;
% plot(traj.z(1,:,cycle_c),traj.z(2,:,cycle_c),'-*b');
% plot(x_pos(1,:,cycle_c), x_pos(2,:,cycle_c),'-*r');
% hold off;
% title('trajectory');
% legend('groundtruth','measurement','GSM');
% xlabel('x');
% ylabel('y');
% 
% z_RMSE = traj.z(:,:,cycle_c)-traj.x(1:2,:,cycle_c);
% z_RMSE = sqrt((z_RMSE(1,:).^2+z_RMSE(2,:).^2)/2);
% Position_RMSE = x_pos(1:2,:,cycle_c)-traj.x(1:2,:,cycle_c);
% Position_RMSE = sqrt((Position_RMSE(1,:).^2+Position_RMSE(2,:).^2)/2);
% subplot(2,2,2);
% plot(log(z_RMSE),'-*b');
% hold on;
% plot(log(real(Position_RMSE)),'-*r');
% hold off;
% title('Error');
% legend('measurement','GSM');
% xlabel('time');
% ylabel('LOG error');
% 
% disp(mean(z_RMSE));
% disp(mean(real(Position_RMSE)));
% 
% subplot(2,2,[3 4]);
% semilogy(traj.s(:,cycle_c),'-og');
% hold on;
% semilogy(1./K_lamda(:,cycle_c),'-or');
% hold off;
% title('univariate part');
% legend('scale','estimation');
% xlabel('time');
% ylabel('value');

end

    if strcmp(state, 'fail')
        RKF   = struct;
        RKF.state   = 'fail';
    else
        RKF   = struct;
        RKF.state   = 'success';
        RKF.x_pos   = x_pos;
    %   RKF.P_pos   = P_pos;
    %   RKF.K_lamda = K_lamda;
        RKF.time    = toc/cycle; % time for every cycle
        RKF.iteration = ite;
        RKF.num_second = num_second;
    end



end

