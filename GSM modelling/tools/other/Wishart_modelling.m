function [dof,V] = Wishart_modelling(x)
% modelling of Wishart modelling function
% dof is the degree of freedom, V is the scale matrix

N = size(x,3);
d = size(x,1);

m_Wishart = mean(x,3); 

f1 = 0;
for i = 1: N
    f1 = f1 + log(det(x(:,:,i)));
end

loss = @(n) -Wishart_loss(x,f1,m_Wishart,n);

% search the dof between 0 and infinity
dof  = fmincon(loss,10,[],[],[],[],d-1,Inf);

V = mean(x,3)/dof; 


% dof_index = 1:1:10;
% for i = 1: size(dof_index,2)
% y(i) = loss(dof_index(i));
% end
% figure(1);
% plot(dof_index,y,'-or');
end

function [out] = Wishart_loss(x,f1,m_Wishart,n)

d = size(x,1);
N = size(x,3);

% mean of wishart distribution

V  = m_Wishart/n;

f2 = 0;
for i = 1: N
    f2 = f2 + trace(V\x(:,:,i));
end

f3 = N*(d*log(2)+log(det(V)));

f4 = 0;
for i = 1: d
    f4 = f4+ gammaln(n/2+(1-i)/2);
end

out = n/2*f1 -0.5*f2 - n/2*f3 - N*f4;
end

