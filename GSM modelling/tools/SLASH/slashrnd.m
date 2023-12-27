function [Y] = slashrnd(mu,sigma,df,cases)
% function for random sampling of slash distribution;
% mu- mean; sigma- correlation matrix; df- dof value;
% cases- number of samples;

d = size(sigma,1);
x = random('Beta',df/2,1,[1,cases]);
Y = zeros(d,cases);
for i =1:cases
    Y(:,i) = mu + (1/sqrt(x(i)).*mvnrnd(zeros(d,1),sigma))';
end
end

