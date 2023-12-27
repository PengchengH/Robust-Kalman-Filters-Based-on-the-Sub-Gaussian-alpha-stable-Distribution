function [pdf] = slashpdf(X, mu, SIGMA, v)
% X samples, mu mean, SIGMA correlation matrix, v dof value;
d  = size(X,1);

R = cholcov(SIGMA);

% disp(size(X));
% disp(mu);
% disp(size(R));
% disp(SIGMA);

f1  = (X-mu)'/R;

f1  = sum(f1.^2,2);

% find zero element and non zero element
ind1 = find(~f1); 
ind2 = find(f1); 

pdf  = zeros(size(X,2),1);
f2   = zeros(size(X,2),1);
f3   = zeros(size(X,2),1);

% for X = mu
pdf(ind1) = v*(det(SIGMA)^(-0.5))/((2*pi)^(d/2))/(v+d);  

% for X != mu
f2(ind2)  = gamma((v+d)/2)*gammainc(f1(ind2)/2,(v+d)/2);
f3(ind2)  = f1(ind2).^((v+d)/2);
pdf(ind2) = v*2^((v+d)/2-1)*det(SIGMA)^(-0.5)/((2*pi)^(d/2)) .*f2(ind2)./f3(ind2);

end

