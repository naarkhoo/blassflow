addpath test_data
load X_nig
K = C*kappa^2 + G;
n  =length(X);
p = -1*ones(n,1);
a = mu^2/sigma^2*ones(n,1) +2;
b = (K*X).^2 + diag(C).*nu;
sim = 1;
[MC_res, numerical_res] = genertate_testtable_rgig(p, a, b, sim);
Q = sum(MC_res(:,1))
Q_num = sum(numerical_res(:,1))
H = (K*X)'*(MC_res(:,2).*(K*X))
H_num = (K*X)'*(numerical_res(:,2).*(K*X))
[EG,EiG, ElogG] = Egig_par(p, a, b);