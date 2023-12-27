classdef vg
   
    properties
        name
        Ke
        fmrdn
        flogpdf
        option
%         frdn
    end

    methods   
       function obj = vg()  
           % stable sub-Gaussian SSG
           obj.name  = 'vg';
           obj.Ke    = @(X_in) X_in;
           obj.fmrdn = @(v_in,s1_in,s2_in) 1./gamrnd(v_in/2,2/v_in,s1_in,s2_in);
           obj.flogpdf  = @(X_in, mu_in, C_in, df_in) vglogpdf(X_in, mu_in, C_in, df_in);
        
       end 
    
    end
    
    methods 
      
        function [out] = opt(~,n,dof,eta1)
            out    = struct;
            out.EK = ( n-dof-2 + sqrt((n-dof-2)^2+4*dof*eta1) )/eta1/2;   
        end           
    
    end
end

