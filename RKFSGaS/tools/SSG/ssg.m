classdef ssg
   
    properties
        name
        Ke
        fmrdn
        gaussleg
        option
    end

    methods   
       function obj = ssg()  
           % stable sub-Gaussian SSG
           obj.name  = 'ssg';
           obj.Ke    = @(X_in) 1./X_in;
           obj.fmrdn = @(alpha_in,s1_in,s2_in) stblrnd(alpha_in/2,1,(cos(pi*alpha_in/4))^(2/alpha_in),0,[s1_in s2_in]);
           load('tools/SSG/gaussleg_para.mat');
           obj.gaussleg= gaussleg;
           
           % Default set
           obj.option.optname      = 'IGGL';
           obj.option.num_series   = 31;
           obj.option.num_point    = 30;
           obj.option.num_particle = 1000;
           
       end 
    
    end
    
    methods 

      function [Y,scale] = frnd(~,mu,sigma,alpha,n)   
          d = size(sigma,1);
          alpha_1 = alpha/2;
          gam     = (cos(pi*alpha/4))^(2/alpha);
          scale       = stblrnd(alpha_1,1,gam,0,[1 n]);
          Y       = zeros(d,n);
          for i =1:n
              Y(:,i) = mu + sqrt(scale(i))*mvnrnd(zeros(d,1),sigma);
          end
      end
      
      function [out] = opt(obj,n,dof,eta1)  
      % The selection of 'optname':
      % -'Gauss Laguerre','Importance Sampling', and 'Inverse Gamma'
      % where 'Gauss Laguerre' is the default strategy
      
      % corresponding property:
      % -'num_point', 'num_particle', 'num_series'
      % default values are: 30; 1000; 31;
      
      % Selecting the specific option follows the examples below.
      % option = struct;
      % option.optname = 'Importance Sampling';
      % option.num_particle = 500;
      % obj.option = option.
 
           opti = obj.option;
           if strcmp(opti.optname,'GL')
               [out] =opt_GL(obj,n,dof,eta1);
           elseif strcmp(opti.optname,'IS')
               [out] =opt_IS(obj,n,dof,eta1);
           elseif strcmp(opti.optname,'IGGL') ||...
                     strcmp(opti.optname,'IGIS')
               [out] =opt_IG(obj,n,dof,eta1); 
           else
               error('option.optname does not exist');
           end  
      end
      
       
      function [out] =opt_GL(obj,n,dof,eta1)
         %default value for the number of series
         num_point = obj.option.num_point;
         address   = find(obj.gaussleg.point_index==num_point);
         
         zerosave = obj.gaussleg.para{1,address}.zero;
         w        = obj.gaussleg.para{1,address}.weight;

         alpha = dof/2;
         beta   = 1;
         gam    = (cos(pi*dof/4))^(2/dof);

         pd1  =  makedist('Stable','alpha',alpha,'beta',1,'gam',gam,'delta',beta*gam*tan(pi*alpha/2));
         pdf1 =  log(pdf(pd1,eta1./zerosave/2));

         f1   =  (n/2-2)*log(zerosave)+pdf1;
         f2   =  log(zerosave) + f1;
         w1   =  log(w);
         
         out = struct;
         out.EK  = sum(exp(w1+f2))/sum(exp(w1+f1))/(eta1/2);

      end
       
      function [out] =opt_IG(obj,n,dof,eta1)
          
          % default value for the number of series
          num_series = obj.option.num_series;

          alpha1 = dof/2;

          C  = zeros(num_series,1);
          m  = zeros(num_series,1);
          cm = zeros(num_series,1);
          b  = eta1/2;

          F_sum = zeros(num_series+1,1);

          % k_inf = min(168,168/alpha1);
          % th = log(2) + log(abs(real(exp(log(gamma(1-k_inf*alpha1))...
          %     +gammaln(k_inf*alpha1+alpha1+n/2+1)-log(k_inf)...
          %     -log(gamma(1-k_inf*alpha1-alpha1))-gammaln(k_inf*alpha1+n/2)))))/alpha1;
          % disp(log(eta1));
          % disp(th);
          % if log(eta1)>th

          for k = 1:num_series-1  
              C(k)  = log(alpha1^2/b) + log((-1)^(k+1)) + gammaln(k*alpha1) -...
                    gammaln(k) -log(gamma(1-k*alpha1)) + log(k/(b^(k*alpha1)));
              C(k)  = exp(C(k));
              m(k)  = (alpha1*k+n/2)/b;
              cm(k) =  C(k)*m(k);
              F_sum(k+1) = F_sum(k) + cm(k);
              if k>4 && sum(abs(real(cm(k-4:k))))/sum(abs(real(F_sum(k-3:k+1))))< 1e-2
                  break;
              end
          end


          if  k < num_series-1
          
              C     = real(C);
              F_sum = real(F_sum);
              EK    = real(F_sum(k+1))/sum(C(1:k+1));
              
              % indicate use the IG strategy
              method_indicator = 0;
          else
              if strcmp(obj.option.optname,'IGGL')
                  out1 =opt_GL(obj,n,dof,eta1);
                  EK   = out1.EK;
              elseif strcmp(obj.option.optname,'IGIS')
                  out1 =opt_IS(obj,n,dof,eta1);  
                  EK   = out1.EK;
              end
              % indicate use the second strategy
              method_indicator = 1;
          end
          
          out = struct;
          out.EK = EK;
          out.indicator = method_indicator;
      end
      
      
      function [out] =opt_IS(obj,n,dof,eta1)
         
         % default value for the number of particles
         num_particle = obj.option.num_particle;

         Kernel = obj.Ke;

         alpha1 = dof/2;
         gam     = (cos(pi*dof/4))^(2/dof);

         x       = stblrnd(alpha1,1,gam,0,[1 num_particle]);
         kx      = Kernel(x);

         weight  = exp(-0.5*eta1.*kx).*(kx.^(n/2));
         weight  = weight/sum(weight);
         
         out = struct;
         out.EK = sum(kx.*weight);

         for i = 1:size(weight,2)
             if weight(i) < 1/ size(x,2)
                 weight(i) = 0;
             end
         end
%          figure(2);
%          plot(x,log(weight),'ob');

       end

         
    
    end
end











