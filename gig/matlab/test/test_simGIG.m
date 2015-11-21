close all
sim = 510000;
p = -1;
%   -0.0007
%    0.4129

b = 10^-1;

a =2;
max_mes = 60;
beta = sqrt(a*b);
alpha = sqrt(a/b)
2/3 * sqrt(1-p)<beta
x = linspace(0,max_mes,1000)';
gam_sim2 = rgig_par(p * ones(sim,1), a* ones(sim,1), b* ones(sim,1),ceil(10^3*rand(6,1)) +1);
f=gig_den(x,p,a,b);
quad(@(x)gig_den(x,p,a,b) ,0,max_mes)
plot(x,f)
 sum( gam_sim2<max_mes)/sim
hold on
plot(x, sqrt(f).*x,'r')
figure(2)
hist(gam_sim2(gam_sim2<max_mes),20)

g = @(x) x.^(p-1) .* exp(-beta/2 * (x + 1./x));
%x0 = 0.0455367
% two_d_beta = 45.7507
% A1 = 0.657937
% A = 7.82264