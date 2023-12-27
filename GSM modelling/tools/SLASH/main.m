clc;
clear;

% X   = random('Beta',1,3,[1,100]);


d     = d;
X     = -20:0.1:20;
mu    = one(d,1);
SIGMA = diag(4*ones(d,1));
v     = 3.5;

[Z] = slashrnd(mu,SIGMA,v,1000);
Y   = slashpdf(X, mu, SIGMA, v);

loss_1 = @(dof_in) -sum(log(slashpdf(Z,mu,SIGMA,dof_in)));
j=0;
XL = v-1;
XU = v+1;
for i = XL:0.1:XU
     j=j+1;
     f(j) = loss_1(i);
end  
figure(1); 
subplot(2,2,1);
semilogy(XL:0.1:XU,f,'-ob','LineWidth',2,'MarkerSize',10);
title('loss of dof');
subplot(2,2,2);
plot(Z,slashpdf(Z, mu, SIGMA, v),'or','LineWidth',2,'MarkerSize',10);
title('loss of dof');

subplot(2,2,3);
plot(X,Y,'Ob');
hold on;
plot(Z,zeros(size(Z)),'Or');
hold off;
