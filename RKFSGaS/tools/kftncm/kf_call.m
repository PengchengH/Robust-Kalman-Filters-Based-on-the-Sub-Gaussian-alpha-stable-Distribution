function [out] = kf_call(traj)

% extract factors
F  = traj.F;
H  = traj.H;
x0 = traj.x0;
P0 = traj.P0;

z  = traj.z;
s  = traj.s;
Qs = traj.Qs;
Rs = traj.Rs;
n  = size(Qs,1);

length = traj.length;
cycle  = traj.cycle;


x_pos = zeros(n,length,cycle);
P_pos = zeros(n,n,length,cycle);

for j = 1:cycle
    x_pos(:,1,j)   = x0(:,j);
    P_pos(:,:,1,j) = P0(:,:,j);
end

state  = 'success';

tic; % timer
for cycle_c = 1:cycle 
    
    for t = 2:length   
        mea   =  z(:,t,cycle_c);
        scale =  s(t,cycle_c);
        xt_   =  x_pos(:,t-1,cycle_c);
        Pt_   =  P_pos(:,:,t-1,cycle_c);
      
        [xt,Pt] = kftncm(xt_,Pt_,mea,scale,Qs,Rs,F,H);
      
        x_pos(:,t,cycle_c)   = xt;
        P_pos(:,:,t,cycle_c) = Pt; 
        
        if sum(isnan(xt),'all')>0 || sum(isnan(Pt),'all')>0
            state = 'fail';
            break;
        end
        
%           disp(Pt);
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
% semilogy(s(:,cycle_c),'-og');
% title('univariate part');
% legend('scale');
% xlabel('time');
% ylabel('value');
    
end
    out   = struct;
    if strcmp(state, 'fail')
        out.state   = state;
    elseif  strcmp(state, 'success')
        out.state   = state;  
        out.x_pos   = x_pos;
        out.time    = toc/cycle;
    else
        error('the state value is illegal'); 
    end

end

