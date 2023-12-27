classdef gm
   
    properties
        options
        Ke
    end

    methods   
    function obj = gm()  
        % default gauss mixture distribution parameter
        obj.options.wp1   = 0.9;
        obj.options.K     = 1;  
        Ke    = @(X_in) 1/X_in;
    end 
    
    end
    
    methods 
    
        function [u,scale] = frnd(obj,Q,U,K,wp)
        % produce two-Gaussian mixture distribution sample
        % Q is the covariance matrix and U*Q is for the larger covariance
        % Q and U are necessary input; K-number of sample, wp1- probablity
        % of the first Gaussian, default value are 1 and 0.9 respectively.
        
            if nargin < 2 || isempty(Q) || isempty(U)
                error(message('gmrnd:TooFewInputs'));
            end

            if isempty(K)
                K = obj.options.K;
            end

            if isempty(wp)
                wp = obj.options.wp1;
            end

            d   = size(Q,1);
            u   = zeros(K,d);
            scale = zeros(K,1);

            for i=1:K
               s_d = rand(1);
               if s_d<wp
                  u(i,:) = mvnrnd(zeros(d,1),Q);
                  scale(i) = 1;
               else
                  u(i,:) = mvnrnd(zeros(d,1),U*Q); 
                  scale(i) = U;
               end 
            end 

            u = u';
        end
    
    end
    
end











