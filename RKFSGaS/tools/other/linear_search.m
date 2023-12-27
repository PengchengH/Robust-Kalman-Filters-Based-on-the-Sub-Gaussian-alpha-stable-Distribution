function [XF] = linear_search(loss,XL,XU,epsilon,iteration)
% here we present an linear search function for one-dimension minimisation.
% loss is the function, VL and VU is the band, epsilon .

for i = 1: iteration
Xa = (XL + XU)/2;
Xb = Xa + epsilon;

Va = loss(Xa);
Vb = loss(Xb);

if Va < Vb
   XU = Xb;
elseif Va > Vb
   XL = Xa;
else
   XU = Xb;  XL = Xa;
end

if XU-XL<= 2*epsilon
    break;
end
end

XF = (XU+XL)/2;

end

