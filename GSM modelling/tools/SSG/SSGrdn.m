function [Y] = SSGrdn(mu,sigma,alpha,n)

d = size(sigma,1);

alpha_1 = alpha/2;
gam     = (cos(pi*alpha/4))^(2/alpha);
x       = stblrnd(alpha_1,1,gam,0,[1 n]);
Y       = zeros(d,n);
for i =1:n
    Y(:,i) = mu + (sqrt(x(i)).*mvnrnd(zeros(d,1),sigma))';
end

end

