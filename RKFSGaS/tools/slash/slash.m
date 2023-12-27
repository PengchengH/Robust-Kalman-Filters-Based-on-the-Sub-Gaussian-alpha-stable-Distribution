classdef slash
   
    properties
        name
        Ke
        fmrdn
        fpdf
        
        option
%         frdn
    end

    methods   
       function obj = slash()  
           obj.name  = 'slash';
           obj.Ke    = @(X_in) X_in;
           obj.fmrdn =  @(v_in,s1_in,s2_in) random('Beta',v_in/2,1,[s1_in,s2_in]);
           obj.fpdf  = @(X_in, mu_in, C_in, df_in) slashpdf(X_in, mu_in, C_in, df_in);
%            obj.frdn  = @(mu_in,C_in,df_in,cases_in) slashrnd(mu_in,C_in,df_in,cases_in);
       end 
    
    end
    
    methods 
        
      function [out] = opt(~,n,dof,eta1)
         % The theoretical soluiton is available
            alpha = (n+dof)/2;
            beta  = eta1/2;
            f1 = (gammainc(beta,alpha+1)*gamma(alpha+1))/(beta^(alpha+1));
            f2 = (gammainc(beta,alpha)*gamma(alpha))/(beta^(alpha));
            out = struct;
            out.EK = f1/f2;
      end 
      
      
      function [out] = opt_IS(obj,n,dof,eta1)
         % default value for the number of particles is 1e3, you can use
         % obj.option.num_particle to adjust the number
          num_particle = obj.option.num_particle; 
 
          Kernel = obj.Ke;
            
          x      = obj.fmrdn(dof,1,num_particle);
          kx     = Kernel(x);

          weight = exp(-0.5*eta1.*kx).*(kx.^(n/2));
          weight = weight/sum(weight);
          
          out    = struct;
          out.EK = sum(kx.*weight);

      end
      
      function [Y,s] = frnd(~,mu,sigma,df,cases)
          % function for random sampling of slash distribution;
          % mu- mean; sigma- correlation matrix; df- dof value;
          % cases- number of samples;
          if (nargin < 3)
               error(message('stats:TooFewInputs'));
          end
          
          if nargin<4 || isempty(cases)
                cases = 1;
          end

          d = size(sigma,1);
          s = random('Beta',df/2,1,[1,cases]);
          Y = zeros(d,cases);
          for i =1:cases
               Y(:,i) = mu + (1/sqrt(s(i)).*mvnrnd(zeros(d,1),sigma))';
          end
      end
    
    end
end

