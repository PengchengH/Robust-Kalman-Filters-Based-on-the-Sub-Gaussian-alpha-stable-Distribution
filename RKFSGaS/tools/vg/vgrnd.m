function [Y] = vgrnd(mu,sigma,df,cases)
% function for random sampling of slash distribution;
% mu- mean; sigma- correlation matrix; df- dof value;
% cases- number of samples;

d = size(sigma,1);

% sampling from gamma and calculate the reciprocal
x = random('Gamma',df/2,2/df,[1,cases]);

x = 1./x;

Y = zeros(d,cases);
for i =1:cases
    Y(:,i) = mu + (1/sqrt(x(i)).*mvnrnd(zeros(d,1),sigma))';
end
end

