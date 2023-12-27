function [out] = assimation(z,xt_,Pt_,para,traj,GSM,varargin)

% default set
ite_num = 50;

if ~isempty(varargin)
   if isstruct(varargin{end})
      option = varargin;
      if ~isempty(option.ite_num)
         ite_num = option.ite_num;         
      end
   else
      error('OPTIONS must be a structure created with STATSET'); 
   end
end

% trajectory Parameter reading
F = traj.F;
H = traj.H;
Qn = traj.Qs;

n = size(Qn,1);
m = size(traj.Rs,1);

% modelling parameter
dof = para.delta;
uk  = para.Wdof;
Uk  = para.Wscale;
% Rn  = Uk/(uk-m-1);

% GsM model functions
Ke  = GSM.Ke;

% Estimate Prior distribution
Pn_pri = F*Pt_*F' + Qn; % norminal prediction covariance

% Initialization
E_lamda      = zeros(1,ite_num);
E_lamda(1)   = 1;
E_K_lamda    = zeros(1,ite_num);
E_K_lamda(1) = Ke(E_lamda(1));  % K(lamda)=1/lamda;

E_R_inv        = zeros(m,m,ite_num);

if uk <= m+1 % give minimial value of uk
    uk = m+1+0.1;
end
E_R_inv(:,:,1) = (uk-m-1)*inv(Uk);

P_pos = zeros(n,n,ite_num);
Rm = zeros(m,m,ite_num);

x_pos = zeros(n,ite_num);
x_pri = zeros(n,ite_num);

num_second = 0;   % number of iterations by second strategy;

for i = 1: ite_num-1
    
    % Calculate q^(i+1)[xk]
    Rm(:,:,i) = inv(E_R_inv(:,:,i))/E_K_lamda(i);
    
    x_pri(:,i+1) = F*xt_;
    K = Pn_pri*H'/(H*Pn_pri*H' + Rm(:,:,i));
    x_pos(:,i+1)   = x_pri(:,i+1) + K*(z-H*x_pri(:,i+1));
    P_pos(:,:,i+1) = (eye(n)-K*H) * Pn_pri;

%     disp(P_pos(:,:,i+1))
    
    % Calculate q^(i+1)[Xi]; q^(i+1)[lamda]
    B = (z-H*x_pos(:,i+1)) * (z-H*x_pos(:,i+1))'+...
        H*P_pos(:,:,i+1)*H';    
    b = z-H*x_pos(:,i+1);
    
    % Based on GH Variance Gamma Distribution lamda;
    eta1 = trace(B*E_R_inv(:,:,i));
    opti  = GSM.opt(m,dof,eta1);
    E_K_lamda(i+1)  =  opti.EK;
    if isfield(opti,'indicator')  
          num_second  = num_second + opti.indicator;
    end
    
    D = E_K_lamda(i+1)*( H*P_pos(:,:,i+1)*H' + b*b' );
    
    % Calculate SIGMA and R;
    ukm = uk + 1;
    Ukm = Uk + D;
    
    E_R_inv(:,:,i+1) = (ukm-m-1)*inv(Ukm);  
    
    
    if i>4 
       
        x_diff = sum(abs(x_pos(:,i-3:i+1)-x_pos(:,i-4:i)),'all'); 
        P_diff = sum(abs(P_pos(:,:,i-3:i+1)-P_pos(:,:,i-4:i)),'all'); 
        K_lamda_diff = sum(abs(E_K_lamda(i-3:i+1)-E_K_lamda(i-4:i)),'all'); 

        x_abs  = sum(abs(x_pos(:,i-3:i+1)),'all');
        P_abs  = sum(abs(P_pos(:,:,i-3:i+1)),'all');
        K_lamda_abs =  sum(abs(E_K_lamda(i-3:i+1)),'all');

        if x_diff/x_abs<1e-2 && P_diff/P_abs<1e-2 && K_lamda_diff/K_lamda_abs<1e-2 
            break;
        end  

    end
    
    
    if sum(isnan(x_pos(:,i+1)),'all')>0 || sum(isnan(P_pos(:,:,i+1)),'all')>0 ...
       ||sum(isnan(E_K_lamda(i+1)),'all')>0 
        break;
    end
    
 
end

out = struct;
out.x_pos   = x_pos(:,i+1);
out.P_pos   = P_pos(:,:,i+1);
out.K_lamda = E_K_lamda(i+1);
out.ite     = i;
out.num_second = num_second;

if sum(isnan(out.x_pos),'all')>0 || sum(isnan(out.P_pos),'all')>0 ...
   ||sum(isnan(out.K_lamda),'all')>0
    out.state = 'fail';
else
    out.state   = 'success';
end

% disp(out.num_second/out.ite);
 
% figure(1);
% subplot(1,2,1);
% plot(x_pos(1,:),'-or');
% subplot(1,2,2);
% plot(E_K_lamda,'-or');

end

