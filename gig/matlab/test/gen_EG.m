%%%
%
% Program to test how well Monte Carlo Expctation (using rgig)
%  compares to regular Expctation for Gig distribution!!


%close all
clearvars
close all
addpath ../
%randseed(1);
%rng(1);
mc_sim  = 10^2;
nu2        = 1;

figure(2)

n = 1000;
p = 2-1*rand(n,1);

b = rand(n,1);

a =rand(n,1);
%lambda_vec = 2-1*rand(n,1);
%psi_vec =rand(n,1);
%chi_vec = rand(n,1);

Eig_sim = zeros(n,1);
 Eg_sim = zeros(n,1);
 Elogg_sim = zeros(n,1);
 Eig_sim2 = zeros(n,1);
 Eg_sim2 = zeros(n,1);
 Elogg_sim2 = zeros(n,1);
 YIGY_sim = 0;
 ELogG_sim = 0;
 EsumIG_sim = 0;
 EsumG_sim = 0;
 YIG_sim    = 0;
for i = 1:mc_sim

  
    gam_sim2 = rgig_par(p, a, b,ceil(10^6*rand(6,1)) );
    Eig_sim2 = Eig_sim2 + 1./gam_sim2;
    Eg_sim2  = Eg_sim2  + gam_sim2; 
    Elogg_sim2 = Elogg_sim2 + log(gam_sim2);
    

end


Eig_sim2 = Eig_sim2 ./mc_sim;
Eg_sim2  = Eg_sim2  ./ mc_sim;
Elogg_sim2 = Elogg_sim2 ./ mc_sim;
tic
[Eig,Eg,Elogg] = Egig(p , a, b);
toc
tic
[Eg2,Eig2,Elogg2] = Egig_par(p , a, b);
toc
max(abs(Eig - Eig2))

subplot(311)
plot(Eig    ,'r.')
hold on
plot(Eig_sim2    ,'blackx')
%plot(1./Cr,'x')
subplot(312)

plot(Eg    ,'r.')
hold on
plot(Eg_sim2    ,'blackx')

subplot(313)
plot(Elogg  ,'r.')
hold on
plot(Elogg_sim2    ,'blackx')
%plot(Cr,'x')
ind_ = sqrt(a.* b) <= min(1/2,2/3 *sqrt(1 -p));
ind_(abs(Eg_sim2-Eg) >2);
