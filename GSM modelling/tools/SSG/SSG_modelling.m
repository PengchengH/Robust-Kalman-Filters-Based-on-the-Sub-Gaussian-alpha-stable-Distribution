classdef SSG_modelling
   
    properties
        model
        options
        initial
    end

    methods   
    function obj = SSG_modelling()  
        
        % search iteration
        obj.options.M   = 200;
        obj.options.M1  = obj.options.M-15;
        
        % number of sampling 
        obj.options.N   = 2000;
        
        % linear search parameter
        obj.options.XL        = 0.01;
        obj.options.XU        = 1.98;
        obj.options.epsilon   = 0.002;
        obj.options.iteration = 10;
        
        % initial value
        obj.initial.mu    = 0;
        obj.initial.scale = 1;
        obj.initial.alpha = 0.1;
        
        
        % stable sub-Gaussian SSG
        obj.model.Ke = @(X_in) 1./X_in;
        obj.model.fmrdn = @(dof_in,s1_in,s2_in) stblrnd(dof_in/2,1,(cos(pi*dof_in/4))^(2/dof_in),0,[s1_in s2_in]);
        obj.model.frdn  = @(mu_in,sigma_in,dof_in,cases_in) dof_SSGrdn(mu_in,sigma_in,dof_in,cases_in);
        obj.model.fitting   = @(Samp,SSG,EM)  SSGfit_EM(Samp,SSG,EM);

    end 
    
    end
    
    methods 
    
    function para = fitting(obj,Samp,model_name)
        if strcmp(model_name,'SSG_MLE')
            para = Mle_fitting(obj,Samp);
        elseif strcmp(model_name,'SSG_EM')
            para = EM_fitting(obj,Samp);
        elseif strcmp(model_name,'SSG_Mellin')
            para = Mellin_fitting(obj,Samp);
        elseif strcmp(model_name,'SSG_ECF')
            para = Ecf_fitting(obj,Samp);
        elseif strcmp(model_name,'SSG_Percentile')
            para = Percentile_fitting(obj,Samp);
        end
    end
    
    function para = Mle_fitting(~,Samp)
        % 	 The FIT function of stable distribution interpolates a precomputed
        %    density table for dof>=0.4 instead of the direct integration method to
        %	 find the maximum likelihood estimation. Estimates from McCulloch's quantile
        %	 method is used as an initial values of the optimization procedure.
        %
        % Reference:
        %   [1] J. P. Nolan (2001) "Maximum likelihood estimation and diagnostics for
        %       stable distributions." Lévy processes. Birkhäuser Boston, p379-400.

        params = stablefit(Samp,0.05); 
        para = struct;
        para.mu      = params(4);
        para.SIGMA   = 2*params(3)^2;
        para.dof   = params(1);
        
        para.dof = max(para.dof,0.02);
        para.dof = min(para.dof,1.98);
        para.SIGMA = max(para.SIGMA,0);
    end
        
    function para = Percentile_fitting(~,Samp)  
        params = stblfit(Samp,'percentile');
        
        %METHOD1 Summary of this method goes here
        para = struct;
        para.mu      = params(4);
        para.SIGMA   = 2*params(3)^2;
        para.dof   = params(1);
        
        para.dof = max(para.dof,0.02);
        para.dof = min(para.dof,1.98);
        para.SIGMA = max(para.SIGMA,0);
    end
        
    function para = Ecf_fitting(~,Samp)  
        optio = statset('MaxIter',100);
        params = stblfit(Samp,'ecf',optio);
        
        %METHOD1 Summary of this method goes here
        para = struct;
        para.mu      = params(4);
        para.SIGMA   = 2*params(3)^2;
        para.dof   = params(1);
        
        para.dof = max(para.dof,0.02);
        para.dof = min(para.dof,1.98);
        para.SIGMA = max(para.SIGMA,0);
    end
        
    function para = Mellin_fitting(obj,Samp)
            % extract sample parameters
            n  = size(Samp,2);
            Y1 = abs(Samp);
            
            logY1 = log(Y1);
            k1    = sum(sum(logY1))/n;
            k2    = sum(sum( (logY1-k1).^2))/n;
            dof = sqrt(2*pi^2/(12*k2-pi^2));
            sigma = exp((dof*k1-(dof-1)*psi(1))/dof);           
            
            %METHOD1 Summary of this method goes here
            para = struct;
            para.mu      = obj.initial.mu;
            para.SIGMA   = 2*sigma^2;
            para.dof   = dof;
            
            para.dof = max(para.dof,0.02);
            para.dof = min(para.dof,1.98);
            para.SIGMA = max(para.SIGMA,0);
            
    end
        
    function para = EM_fitting(obj,Samp)     
         % number of iterations
         M  = obj.options.M;
         M1 = obj.options.M1;
         M2 = M-M1;
         N  = obj.options.N;
          
         % extract model parameters
         Ke    = obj.model.Ke;
         fmrdn = obj.model.fmrdn;
        
         % extract sample parameters
         n = size(Samp,2);
         d = size(Samp,1);
         Y = Samp;

         % give initial parameter value by Mellin transform
         mu        = zeros(d,M); % mean
         SIGMA     = zeros(d,d,M); %covariance
         dof        = zeros(1,M); % dof
         
%          XL = obj.options.XL;
%          XU = obj.options.XU;
%          epsilon = obj.options.epsilon;
%         iteration = obj.options.iteration;
         
         ini   = Mellin_fitting(obj,Samp);
         mu(:,1)      = ini.mu;
         SIGMA(:,:,1) = ini.SIGMA;
         dof(1)       = ini.dof; 
                 
%          D = zeros(1,n);
%          et1 = zeros(1,n);
         
         for t = 1:M-1 
            %% Sampling univariate dof stable
             LAMBDA  = fmrdn(dof(t),1,N);
                

            %%  Update the mean estimation;
            %   Calculate Di
%                 for i = 1:n  
%                     D(i) =(Y(:,i)-mu(:,t))'/SIGMA(:,:,t)*(Y(:,i)-mu(:,t));
%                 end
%             
%                 for i = 1:n
%                    et1(i) = sum(Ke(LAMBDA).^(d/2+1).*exp(-D(i)/2.*Ke(LAMBDA)))/...
%                              sum(Ke(LAMBDA).^(d/2).*exp(-D(i)/2.*Ke(LAMBDA)));
%                 end
%                 
%                 [~, col] = find(isnan(et1));
%                 et1(col) = 0;  
%                 
%                 mu(:,t+1) = sum(Y.*et1,2)/sum(et1);
              mu(:,t+1) = zeros(d,1);

            %% Sampling Wi;
                Y1 = Y- mu(:,t+1);
                Y2 = Y1./sqrt(exprnd(1,1,n));

                w  = zeros(1,n);

                for i = 1:n   
                       j = 0;
                       b = d^(d/2)*(Y2(:,i)'/SIGMA(:,:,t)*Y2(:,i))^(-d/2)*exp(-d/2)/...
                         ((2*pi)^(d/2)*det(SIGMA(:,:,t))^0.5);
                       k = 0;
                      while j==0  
%                            D  = Y2(:,i)'/SIGMA(:,:,t)*Y2(:,i);
%                            if (D-d)/d >1e20
%                                j = 1;
%                                w(i) = sqrt(d)*(Y2(:,i)'/SIGMA(:,:,t)*Y2(:,i))^(-0.5);
%                            end
                          
                           w(i) = wblrnd(1,dof(t));        
                           u  = b*rand;         
                           th = w(i)^d*exp(-0.5*(Y2(:,i)'/SIGMA(:,:,t)*Y2(:,i))*w(i)^2)/...
                                ((2*pi)^(d/2)*det(SIGMA(:,:,t))^0.5);
                           if u < th
                               j = 1;
                           end
                           
                           k = k+1;  
                           if k>1e5
                               w(i) = sqrt(d)*(Y2(:,i)'/SIGMA(:,:,t)*Y2(:,i))^(-0.5);
                               j = 1;
                           end 
                               
                      end 
                end

                for i = 1: n
                    SIGMA(:,:,t+1) = SIGMA(:,:,t+1)+ Y2(:,i)*Y2(:,i)'*w(i)^2;
                end
                SIGMA(:,:,t+1) = SIGMA(:,:,t+1)/n;

                cost = @(dof_in) abs(n/dof_in+sum(log(w))-sum(w.^dof_in.*log(w)));
                
                lb = 0.01;
                ub = 1.98;
                dof(t+1)    = fmincon(cost,dof(t),[],[],[],[],lb,ub);
%               dof(t+1) = linear_search(cost,XL,XU,epsilon,iteration);
                
        %% Convergence analysis
        %  &&sum(abs(mu(:,t-M2+1:t+1)-mu(:,t-M2:t)),'all')/sum(abs(mu(:,t-M2+1:t+1)),'all')< 1e-2...
        if t>M2 && sum(abs(dof(t-M2+1:t+1)-dof(t-M2:t)))/sum(dof(t-M2+1:t+1))< 3e-2...
           && sum(abs(SIGMA(:,:,t-M2+1:t+1)-SIGMA(:,:,t-M2:t)),'all')/sum(abs(SIGMA(:,:,t-M2+1:t+1)),'all')< 3e-2
             break;
        end

%             % Display
%                 j=0;
%                 for i = 0.2:0.1:1.9
%                     j=j+1;
%                     f(j) = cost(i);
%                 end
%                 figure(1); 
%                 subplot(2,2,1);
%                 plot(0.2:0.1:1.9,f,'-ob','LineWidth',2,'MarkerSize',10);
%                 title('loss');   
%                 subplot(2,2,2);
%                 plot(mu(1:1:t+1),'-om','LineWidth',2,'MarkerSize',3);
%                 title('mean');
%                 subplot(2,2,3);
%                 plot(squeeze(SIGMA(1,1,1:t+1)),'-om','LineWidth',2,'MarkerSize',3);
%                 title('variance');
%                 subplot(2,2,4);
%                 plot(dof(1:t+1),'-om','LineWidth',2,'MarkerSize',3);
%                 title('dof');
%                 pause(0.1);

        end

        para = struct;
        para.mu    = mean(mu(:,t-M2+1:t+1),2);
        para.SIGMA = mean(SIGMA(:,:,t-M2+1:t+1),3);
        para.dof   = mean(dof(t-M2+1:t+1));
        
        para.dof = max(para.dof,0.02);
        para.dof = min(para.dof,1.98);
        para.SIGMA = max(para.SIGMA,0);
    end 
    
    end
    
    %===========Stable Fitting Functions===========
end

function [parmhat] = stablefit(x,dof)
%STABLEFIT Parameter estimates and confidence intervals for stable data.
% This method is only suitable for the shape parameter dof >= 0.4.

% Reference:
%   [1] J. P. Nolan (2001) "Maximum likelihood estimation and diagnostics for
%       stable distributions." Lévy processes. Birkhäuser Boston, p379-400.

if ~isvector(x)
    error(message('stats:addstable:VectorRequired'));
end

if nargin < 2 || isempty(dof)
    dof = 0.05;
end

% Compute initial estimates by McCulloch quantile method to speed up the
% optimization procedure.
params0    = intMle(x);

% Based on stable sub-Gaussian, we directly give beta = 0; delat =0; 
params0(2) = 0;  
params0(4) = 0;

% If the initial estimates of dof or BETA are near or at the boundary,
% fminsearch may not make any update on the parameter estimates. So, try to
% reassign some reasonable closer values to the initial estimates of dof
% and BETA.
if (2-params0(1) < 1e-3)
    params0(1) = 1.95;
end
% if (1-abs(params0(2)) < 1e-3)
%     params0(2)= sign(params0(2))*0.95;
% end

% Transform the parameters to eliminate bounds in order to apply fminsearch
phi0 = varTrans(params0,'forward');

% Maximize the log-likelihood with respect to the transformed parameters
% When possible, use 'fmincon' for speed
lb = [0.5 0.1];
ub = [1.9  100];

% options = optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon);
% parmhat = particleswarm(@(param) stable_nloglf(x,[param(1),phi0(2),param(2),phi0(4)]),...
%     2 ,lb,ub,options);
[parmhat] = fmincon(...
    @(param)stable_nloglf(x,[param(1),phi0(2),param(2),phi0(4)]),...
     [phi0(1),phi0(3)],[],[],[],[],lb,ub);


parmhat(1) = parmhat(1);
parmhat(3) = parmhat(2);
parmhat(2) = 0;
parmhat(4) = 0;


end

function params = intMle(x)
%INTMLE Parameter estimation of stable data by McCulloch's quantile method.
% This method is only suitable for dof>0.5.

% Reference:
%   [1] J. H. McCulloch (1986) "Simple Consistent Estimators of Stable
%       Distribution Parameters". Cummun. Stat. -Simula., 15(4).

if ~isvector(x)
    error(message('stats:addstable:VectorRequired'));
end

Xpcts = prctile(x,[95,75,50,25,5]);
nudof = min(25,max((Xpcts(1) - Xpcts(5))/(Xpcts(2) - Xpcts(4)),2.439));
nuBeta = min(1,max((Xpcts(1) + Xpcts(5) - 2*Xpcts(3))/(Xpcts(1) - Xpcts(5)),-1));

% Input tables for dof and Beta estimation
nuA = [2.439 2.5 2.6 2.7 2.8 3.0 3.2 3.5 4.0 5.0 6.0 8.0 10 15 25];
nuB = [0 0.1 0.2 0.3 0.5 0.7 1];
[na, nb] = meshgrid( nuA , nuB );
dofTbl=  [2.000 2.000 2.000 2.000 2.000 2.000 2.000;...
    1.916 1.924 1.924 1.924 1.924 1.924 1.924;...
    1.808 1.813 1.829 1.829 1.829 1.829 1.829;...
    1.729 1.730 1.737 1.745 1.745 1.745 1.745;...
    1.664 1.663 1.663 1.668 1.676 1.676 1.676;...
    1.563 1.560 1.553 1.548 1.547 1.547 1.547;...
    1.484 1.480 1.471 1.460 1.448 1.438 1.438;...
    1.391 1.386 1.378 1.364 1.337 1.318 1.318;...
    1.279 1.273 1.266 1.250 1.210 1.184 1.150;...
    1.128 1.121 1.114 1.101 1.067 1.027 0.973;...
    1.029 1.021 1.014 1.004 0.974 0.935 0.874;...
    0.896 0.892 0.887 0.883 0.855 0.823 0.769;...
    0.818 0.812 0.806 0.801 0.780 0.756 0.691;...
    0.698 0.695 0.692 0.689 0.676 0.656 0.595;...
    0.593 0.590 0.588 0.586 0.579 0.563 0.513]';
betaTbl=  [ 0.000 2.160 1.000 1.000 1.000 1.000 1.000;...
    0.000 1.592 3.390 1.000 1.000 1.000 1.000;...
    0.000 0.759 1.800 1.000 1.000 1.000 1.000;...
    0.000 0.482 1.048 1.694 1.000 1.000 1.000;...
    0.000 0.360 0.760 1.232 2.229 1.000 1.000;...
    0.000 0.253 0.518 0.823 1.575 1.000 1.000;...
    0.000 0.203 0.410 0.632 1.244 1.906 1.000;...
    0.000 0.165 0.332 0.499 0.943 1.560 1.000;...
    0.000 0.136 0.271 0.404 0.689 1.230 2.195;...
    0.000 0.109 0.216 0.323 0.539 0.827 1.917;...
    0.000 0.096 0.190 0.284 0.472 0.693 1.759;...
    0.000 0.082 0.163 0.243 0.412 0.601 1.596;...
    0.000 0.074 0.147 0.220 0.377 0.546 1.482;...
    0.000 0.064 0.128 0.191 0.330 0.478 1.362;...
    0.000 0.056 0.112 0.167 0.285 0.428 1.274]';

dof = interp2(na,nb,dofTbl,nudof,abs(nuBeta));

Beta = sign(nuBeta) * interp2(na,nb,betaTbl,nudof,abs(nuBeta));

% Reset dof if necessary.
if dof > 2
    dof = 2;
    Beta = sign(nuBeta);
end

% Reset Beta if necessary.
if Beta > 1
    Beta = 1;
elseif Beta < -1
    Beta = -1;
end

% Input tables for Gam and Delta estimation
va = [2 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5];
vb = [0 0.25 0.5 0.75 1];
gamTbl = [  1.908 1.908 1.908 1.908 1.908;...
    1.914 1.915 1.916 1.918 1.921;...
    1.921 1.922 1.927 1.936 1.947;...
    1.927 1.930 1.943 1.961 1.987;...
    1.933 1.940 1.962 1.997 2.043;...
    1.939 1.952 1.988 2.045 2.116;...
    1.946 1.967 2.022 2.106 2.211;...
    1.955 1.984 2.067 2.188 2.333;...
    1.965 2.007 2.125 2.294 2.491;...
    1.980 2.040 2.205 2.435 2.696;...
    2.000 2.085 2.311 2.624 2.973;...
    2.040 2.149 2.461 2.886 3.356;...
    2.098 2.244 2.676 3.265 3.912;...
    2.189 2.392 3.004 3.844 4.775;...
    2.337 2.635 3.542 4.808 6.247;...
    2.588 3.037 4.534 6.636 9.144]';
deltaTbl = [0.000  0.000  0.000  0.000  0.000;...
    0.000 -0.017 -0.032 -0.049 -0.064;...
    0.000 -0.030 -0.061 -0.092 -0.123;...
    0.000 -0.043 -0.088 -0.132 -0.179;...
    0.000 -0.056 -0.111 -0.170 -0.232;...
    0.000 -0.066 -0.134 -0.206 -0.283;...
    0.000 -0.075 -0.154 -0.241 -0.335;...
    0.000 -0.084 -0.173 -0.276 -0.390;...
    0.000 -0.090 -0.192 -0.310 -0.447;...
    0.000 -0.095 -0.208 -0.346 -0.508;...
    0.000 -0.098 -0.223 -0.383 -0.576;...
    0.000 -0.099 -0.237 -0.424 -0.652;...
    0.000 -0.096 -0.250 -0.469 -0.742;...
    0.000 -0.089 -0.262 -0.520 -0.853;...
    0.000 -0.078 -0.272 -0.581 -0.997;...
    0.000 -0.061 -0.279 -0.659 -1.198]';
[ma, mb] = meshgrid( va, vb );
nuGam = interp2(ma,mb,gamTbl,dof,abs(Beta));
Gam = (Xpcts(2)-Xpcts(4))/nuGam;

% Reset the value of Gam if necessary
if Gam==0 || isnan(Gam) || isinf(Gam)
    s = std(x);
    if s>0
        Gam = s;
    else
        Gam = eps;
    end
end

nuKsi = sign(Beta) * interp2(ma,mb,deltaTbl,dof,abs(Beta));
Ksi = Xpcts(3) + Gam*nuKsi;

if dof == 1
    Delta = Ksi;
else
    Delta = Ksi - Beta*Gam*tan(pi*dof/2);
end
params(1) = dof;
params(2) = Beta;
params(3) = Gam;
params(4) = Delta;
end

function y = tabpdf(x,dof,beta)
%TABPDF Density for standardized stable data by interpolating tabulated density table.
x = atan(x);

% Use cubic interpolation and enforce no extrapolation (this makes
% extrapolation NaN's)
% Keep the griddedInterpolant 'persistent' for speed with repeated function
% calls
persistent G;
if (isempty(G))
    s = load('private/StablePdfTable.mat');
    G = griddedInterpolant({s.b, s.a, s.xgd}, s.p, 'linear','none');
end

y = G({beta,dof,x});
y = reshape(y,size(x));

% Assign 0 to NaN's caused by extrapolation.
y(isnan(y)) = 0;

end

function nll = neglog_pdf(x,dof,beta)
%NEGLOG_PDF Negative log-likelihood for standardized data using interpolated densities.

nll = -log(max(realmin,tabpdf(x,dof,beta))); % In case of small or negative densities
end


function nll = stable_nloglf(x,params)
%STABLE_NLOGLF Objective function for stable maximum likelihood.
n = numel(x);
% params = varTrans(params,'backward');
dof = params(1);
beta = params(2);
gam = params(3);
delta = params(4);
z = (x-delta)./gam;
nll = sum(neglog_pdf(z,dof,beta)) + n*log(gam);
end

function phi = varTrans(theta,direct)
%VARTRANS Transform the variables to apply the optimization of the log-likelihood function.
%
% When DIRECTION = 'FORWARD', THETA is a vector of parameter
% values in the stable distribution and PHI is a vector of transformed
% parameter values.
% When DIRECTION = 'BACKWARD', THETA is a vector of transformed parameter
% and PHI is a vector of the parameter values used in stable distribution.
%
% The optimization algorithm fminsearch is for unconstrained variables. The
% variables THETA=(dof,BETA,GAM,DELTA) in stable distribution are constrained
% such that 0<dof<=2, -1<=BETA<=1, GAM>0 and -Inf<DELTA<Inf. To eliminate the
% bounds of THETA, we can do the following transformation:
%       dof = a + (b - a)./(1 + exp(-PHI(1))), where a=0, b=2
%       BETA = a + (b - a)./(1 + exp(-PHI(2))), where a=-1, b=1
%       GAM = exp(PHI(3))
%       DELTA = PHI(4)
% then the parameters PHI are unbounded.

narginchk(2,Inf);
phi = ones(size(theta));

if strcmpi(direct,'forward')
    phi(1) = -log(2./theta(1) - 1);
    phi(2) = -log((2./(theta(2)+1) - 1));
    phi(3) = log(theta(3));
    phi(4) = theta(4);
else % direct = backward
    phi(1) = 2./(1 + exp(-theta(1)));
    phi(2) = -1 + 2./(1 + exp(-theta(2)));
    phi(3) = exp(theta(3));
    phi(4) = theta(4);
end
end
      









