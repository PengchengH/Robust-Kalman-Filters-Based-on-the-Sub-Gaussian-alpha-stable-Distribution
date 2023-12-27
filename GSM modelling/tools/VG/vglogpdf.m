function [pdf] = vglogpdf(X, mu, SIGMA, v)
% X samples, mu mean, SIGMA correlation matrix, v dof value;
d = size(X,1);
R = cholcov(SIGMA);

f1 = (X-mu)'/R;
f1 = v*sum(f1.^2,2);

% f2 = 2^(1-v/2-d/2)*v^(d/2)*det(SIGMA)^(-0.5)*pi^(-d/2)/gamma(v/2);
% pdf = f2 * besselk(v/2-d/2,sqrt(f1)).* (f1.^((v-d)/4));

f2 = log(2^(1-v/2-d/2)*v^(d/2)*det(SIGMA)^(-0.5)*pi^(-d/2))-gammaln(v/2);
pdf = f2 + log(besselk(v/2-d/2,sqrt(f1),1)) -sqrt(f1) + log(f1.^((v-d)/4));

% pdf = exp(pdf);

end

