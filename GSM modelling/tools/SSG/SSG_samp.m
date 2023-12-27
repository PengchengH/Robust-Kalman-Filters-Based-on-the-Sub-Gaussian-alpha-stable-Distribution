function [Y] = SSG_samp(MD)
    
     sigma = MD.Q;
     alpha = MD.alpha;
     n     = MD.n;
     d     = size(sigma,1);
     mu    = zeros(d,1);

    [Y] = SSGrdn(mu,sigma,alpha,n);
                    
end

