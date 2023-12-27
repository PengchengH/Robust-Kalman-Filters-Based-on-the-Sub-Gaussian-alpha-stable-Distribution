function [model_traj] = tracking_model(mea_noise,dof,model_traj,cycle,varargin)

if strcmp(mea_noise,'gm') && isstruct(varargin{end})
     opt = varargin{end};
     wp  = opt.wp;
elseif strcmp(mea_noise,'gm') && ~isstruct(varargin{end}) 
     error(message(' TooFewInputs'));
end

% extract parameters
length = model_traj.length;
Qs     = model_traj.Qs;
Rs     = model_traj.Rs;
for_op = model_traj.for_op;
obs_op = model_traj.obs_op;

n      = size(Qs,1);
m      = size(Rs,1);

% produce signal and measurement noise and scale ground truth
% for measurement noise

% signal noise
w = zeros(size(Qs,1),length,cycle);
for cycle_c = 1:cycle
    w(:,:,cycle_c) = mvnrnd(zeros(n,1),Qs,length)';
end

% measurement noise
v = zeros(size(Rs,1),length,cycle);
% scale value, scale*Qn is the real Q
s = zeros(length,cycle);
if strcmp(mea_noise,'gm')
    Sampler = gm;
    for cycle_c = 1:cycle
        [out1,out2] = Sampler.frnd(Rs,dof,length,wp);
        v(:,:,cycle_c) = out1;
        s(:,cycle_c)   = out2;
    end
    
elseif strcmp(mea_noise,'stdt')
    Sampler = stdt;
    for cycle_c = 1:cycle
        [out1,out2] = Sampler.frnd(zeros(1,m),Rs,dof,length);  
        v(:,:,cycle_c) = out1';
        s(:,cycle_c)   = out2;
    end
elseif strcmp(mea_noise,'ssg')
    Sampler = ssg;
    for cycle_c = 1:cycle
        [out1,out2] = Sampler.frnd(zeros(1,m),Rs,dof,length);  
        v(:,:,cycle_c) = out1;
        s(:,cycle_c)   = out2;
    end 
end

% Ground Truth Trajectory
x = zeros(size(Qs,1),length,cycle);
z = zeros(size(Rs,1),length,cycle);
for j = 1:cycle
    
    x(:,1,j) = [0,0,10,10]';  
    for i = 2:length
        x(:,i,j) = for_op(x(:,i-1,j),w(:,i,j));
    end  
    z(:,:,j) = obs_op(x(:,:,j),v(:,:,j));
end

% Initial information
P0 = zeros(size(Qs,1),size(Qs,1),cycle);
x0 = zeros(size(Qs,1),cycle);

for j = 1:cycle
   P0(:,:,j) = diag([25 25 2 2]); % intial estimation covariance
   x0(:,j)   = mvnrnd(x(:,1,j),P0(:,:,j)); % initial position
end


model_traj.length = length;
model_traj.x   = x;
model_traj.z   = z;
model_traj.s   = s;
model_traj.P0  = P0;
model_traj.x0  = x0;
model_traj.cycle = cycle;

% for j = 1:cycle
%     figure(1);
%     subplot(2,3,1);
%     plot(w(1,:,j),w(2,:,j),'Or');
%     title('process noise');
%     subplot(2,3,2);
%     loglog(abs(v(1,:,j)),abs(v(2,:,j)),'Or');
%     title('log(abs(measurement noise))');
%     subplot(2,3,3);
%     semilogy(s(:,j),'Ob');
%     title('log(Scale)');
%     subplot(2,3,4);
%     plot(x(1,:,j),x(2,:,j),'-Ob');
%     title('Trajectory');
%     subplot(2,3,5);
%     plot(z(2,:,j),'-Ob');
%     title('Measurement-Y');
%     subplot(2,3,6);
%     plot(x(3,:,j),x(4,:,j),'-Ob');
%     title('Velocity');

%     figure(1);
%     plot(x(1,1:1:50,j),x(2,1:1:50,j),'.b','LineWidth',2,'MarkerSize',10);
%     hold on;
%     plot(z(1,1:1:50,j),z(2,1:1:50,j),'.r','MarkerSize',10);
%     legend('X_k','Z_k');
%     xlabel('position-x(m)');
%     ylabel('position-y(m)');

% end


end



