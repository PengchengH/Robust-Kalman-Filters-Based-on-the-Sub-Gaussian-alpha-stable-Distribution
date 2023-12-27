% Tracking model initial parameters
T = 1; % sampling interval
length = 300;

% state transition and measurement matrix
F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1]; 
H = [1 0 0 0; 0 1 0 0]; 

% signal and measurement noise parameter
Qs = 0.1*[T^3/3 0 T^2/2 0; 0 T^3/3 0 T^2/2;...
          T^2/2 0 T     0; 0 T^2/2 0 T];
Rs = ratio_scale*10*eye(2);
 
% forcast and observation
for_op = @(x,w) F*x+w;
obs_op = @(x,v) H*x+v;

traj   = struct;
traj.T = T; % sampling interval
traj.F = F;
traj.H = H;

traj.length = length;

traj.Qs = Qs;
traj.Rs = Rs;

traj.for_op = for_op;
traj.obs_op = obs_op;

