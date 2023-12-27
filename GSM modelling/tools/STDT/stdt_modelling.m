classdef stdt_modelling
   
    properties
        model
        options
        initial
    end

    methods   
    function obj = stdt_modelling()  
        
        % search iteration
        obj.options.M   = 300;
        obj.options.M1  = obj.options.M-5;
        
        % number of sampling 
        obj.options.N   = 2000;
        
        % linear search parameter
        obj.options.XL        = 0.01;
        obj.options.XU        = 20;
        obj.options.epsilon   = 0.01;
        obj.options.iteration = 12;
        
        % initial value
        obj.initial.mu    = 0;
        obj.initial.scale = 1;
        obj.initial.dof   = 1;
        
        
        % stable sub-Gaussian SSG
        obj.model.Ke    = @(X_in) X_in;
        obj.model.fmrdn = @(v_in,s1_in,s2_in) gamrnd(v_in/2,2,s1_in,s2_in)./v_in;
        obj.model.fpdf  = @(X_in, mu_in, C_in, df_in) stpdf(X_in, mu_in, C_in, df_in);
        obj.model.frdn   = @(mu_in,C_in,df_in,cases_in) strnd(mu_in',C_in,df_in,cases_in)';

    end 
    
    end
    
    methods 
    
    function para = fitting(obj,Samp,model_name)
        if strcmp(model_name,'stdt_EM')
            para = EM_fitting(obj,Samp);
        elseif strcmp(model_name,'stdt_MLE')
            para = MLE_fitting(obj,Samp);
        end
    end
    
    function [para] = MLE_fitting(obj,Samp)
        fpdf  = obj.model.fpdf;
        d     = size(Samp,1);
        Y     = Samp;
        if  d>1
        error('dimension shoule be 1');
        end
        
        mu     = obj.initial.mu*ones(d,1);
        SIGMA  = obj.initial.scale*diag(ones(d,1));
        dof      = obj.initial.dof;
        
%       loss = @(dof_in,SIGMA) -sum(log(fpdf(Y,zeros(d,1),SIGMA,dof_in)));
        loss = @(x) -sum(log(fpdf(Y,mu,x(1),x(2))));
        x0   = [SIGMA,dof];
        lb    = [0 0];
        ub    = [Inf  Inf];
        
%       options = optimset('PlotFcns',@optimplotfval);
        x    = fmincon(loss,x0,[],[],[],[],lb,ub);

        para = struct;
        para.mu    = mu;
        para.SIGMA = x(1);
        para.dof   = x(2);

    end
    
    function [para] = EM_fitting(obj,Samp)
        % number of iterations
        M  = obj.options.M;
        M1 = obj.options.M1;
        M2 = M - M1;
        N  = obj.options.N;
          
        % extract model parameters
        Ke    = obj.model.Ke;
        fmrdn = obj.model.fmrdn;
        fpdf  = obj.model.fpdf;
        
        % extract sample parameters
        n = size(Samp,2);
        d = size(Samp,1);
        Y = Samp;
        
        % give initial parameter value.
        mu           = zeros(d,M); % mean
        SIGMA        = zeros(d,d,M); %covariance
        dof          = zeros(1,M); % dof
        mu(:,1)      = obj.initial.mu*ones(d,1);
        SIGMA(:,:,1) = obj.initial.scale*diag(ones(d,1));
        dof(1)       = obj.initial.dof;
        
        % linear search scale
        XL        = obj.options.XL;
        XU        = obj.options.XU;
        epsilon   = obj.options.epsilon;
        iteration = obj.options.iteration;

        D      = zeros(1,n);
        et1    = zeros(1,n);

        for t = 1:M-1 
        % Sampling univariate mixture distribution
        LAMBDA  = fmrdn(dof(t),1,N);

        %%  Update the mean estimation;
        %   Calculate Di
%         for i = 1:n  
%         D(i) =(Y(:,i)-mu(:,t))'*inv(SIGMA(:,:,t))*(Y(:,i)-mu(:,t));
%         end
% 
%         for i = 1:n
%         et1(i) = sum(Ke(LAMBDA).^(d/2+1).*exp(-D(i)/2.*Ke(LAMBDA)))/...
%              sum(Ke(LAMBDA).^(d/2).*exp(-D(i)/2.*Ke(LAMBDA)));
%         end
% 
%         [~, col] = find(isnan(et1));
%         et1(col) = 0;   
% 
%         mu(:,t+1) = sum(Y.*et1,2)/sum(et1);
        mu(:,t+1) = zeros(d,1);

        %%  Update the SIGMA estimation;    
        for i = 1:n  
        D(i) =(Y(:,i)-mu(:,t+1))'*inv(SIGMA(:,:,t))*(Y(:,i)-mu(:,t+1));
        end

        for i = 1:n
        et1(i) = sum(Ke(LAMBDA).^(d/2+1).*exp(-D(i)/2.*Ke(LAMBDA)))/...
             sum(Ke(LAMBDA).^(d/2).*exp(-D(i)/2.*Ke(LAMBDA)));
        end

        [~, col] = find(isnan(et1));
        et1(col) = 0;  

        for i = 1:n
        SIGMA(:,:,t+1) = SIGMA(:,:,t+1)+(Y(:,i)-mu(:,t+1))*(Y(:,i)-mu(:,t+1))'*et1(i);
        end
        SIGMA(:,:,t+1) = SIGMA(:,:,t+1) /n;

        %%  Update the def estimation;                         
        loss_1 = @(dof_in) -sum(log(fpdf(Y',mu(:,t+1)',SIGMA(:,:,t+1),dof_in)));
        dof(t+1) = linear_search(loss_1,XL,XU,epsilon,iteration);
        
        %% Convergence analysis
        %  &&sum(abs(mu(:,t-M2+1:t+1)-mu(:,t-M2:t)),'all')/sum(abs(mu(:,t-M2+1:t+1)),'all')< 1e-2...
        if t>M2 && sum(abs(dof(t-M2+1:t+1)-dof(t-M2:t)))/sum(dof(t-M2+1:t+1))< 1e-2...
           && sum(abs(SIGMA(:,:,t-M2+1:t+1)-SIGMA(:,:,t-M2:t)),'all')/sum(abs(SIGMA(:,:,t-M2+1:t+1)),'all')< 1e-2
             break;
        end

        %%  Display
            j=0;
            for i = XL:0.3:XU
                j=j+1;
                f(j) = loss_1(i);
            end  
            figure(1); 
            subplot(2,2,1);
            plot(XL:0.3:XU,f,'-ob','LineWidth',2,'MarkerSize',10);
            title('loss of dof');
            subplot(2,2,2);
            plot(mu(1:t+1),'-om','LineWidth',2,'MarkerSize',10);
            title('mean');
            subplot(2,2,3);
            plot(squeeze(SIGMA(1,1,1:t+1)),'-om','LineWidth',2,'MarkerSize',10);
            title('variance');
            subplot(2,2,4);
            plot(dof(1:t+1),'-om','LineWidth',2,'MarkerSize',10);
            title('degree of freedom');
            pause(0.1);

        end

        para = struct;
        para.mu    = mean(mu(:,t-M2+1:t+1),2);
        para.SIGMA = mean(SIGMA(:,:,t-M2+1:t+1),3);
        para.dof   = mean(dof(t-M2+1:t+1));

    end

    
    end
end











