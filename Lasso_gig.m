% baysian lasso
% based on Korobilis, D  Hierarchical Shrinkage Priors for Dynamic Regressions with Many Predictors. International Journal of Forecasting 2013 implementation
% https://sites.google.com/site/dimitriskorobilis/matlab/bayesian-lasso

clear all; clc;

N = 100;
p = 128;

load('simmat.mat')

D = simmat;
D_normalized = (D/max(D(:))).^2; % should it be normalize like this ?

x=randn(N,p)*chol(abs(1-D_normalized));
beta_true = [zeros(1, 64), ones(1, 64)];
sigma_sim = 1;
y = x * beta_true' + 1*sigma_sim*randn(N,1);

% Standardize the data and demean y as it is required
x = standardize(x);
y = y-mean(y);

% ----------------Gibbs related preliminaries
nsave = 1000;
nburn = 1000;
ntot = nsave + nburn;

beta_draws = zeros(nsave,p);
tau2_draws = zeros(nsave,p);
sigma2_draws = zeros(nsave,1);
% ----------------Set priors

tau_prior_gamma_c = exprnd(1); %1 %page 179, chapter 3.3
if p > (N-1) % page 179, chapter 3.3
    beta_ls = lscov(x, y);
    M = 1/N * sum(beta_ls.^2);  
    beta = 0*ones(p,1);
else
    M = 1/p * sum(beta_OLS.^2);
    beta = inv(x'*x)'*(x'*y);
end
tau_prior_gamma_d = sqrt(Draw_iGamma(2, M)/(2*tau_prior_gamma_c));  %1
% http://emmpa2014.rimir.ro/attachments/article/20/keynote2-slides.pdf
% suggest both parameters less than 1

% Initialize parameters
tau2 = 4*ones(p,1);
V_L = diag(tau2);
lambda2 = 0.1; %(p*sqrt(sigma2_OLS)/(sum(abs(beta_OLS)))).^2;
sigma2 = 0.1; %sigma2_OLS;


%==========================================================================
%====================| GIBBS ITERATIONS START HERE |=======================
disp('Now you are running the Bayesian Lasso with IG prior')
if_accepted = zeros(1, ntot);
acceptance_probability_lambda_gamma_vec = zeros(1, ntot);
for irep = 1:ntot
    irep
    if mod(irep,500)==0
        disp(['Iteration: ' num2str(irep) ' of ' num2str(ntot) '. Please wait a few seconds...'])
    end
    
    beta_old = beta;
    
    % 1. Update beta from Normal
    A = inv(x'*x + inv(V_L ));
    post_mean_beta = A*x'*y; 
    post_var_beta = sigma2*A;
    beta = mvnrnd(post_mean_beta, post_var_beta)';
    
    % 2. Update tau2_j from GIG
    for j = 1:p
        tau_post_p = tau_prior_gamma_c - 0.5;
        tau_post_b = beta(j,1).^2 / sigma2; % b/x
        tau_post_a = tau_prior_gamma_d * 2; % ax
        tau2_inverse = gig_matlab(tau_post_p, tau_post_a, tau_post_b);
        tau2(j,1) = 1/tau2_inverse;
    end

    tune_gam = 1.2; %20-30 acceptance
    z = randn(1);
    gam_lam = tau_prior_gamma_c;
    gam_lam_proposal = exp(tune_gam * z) * gam_lam;
    
    p_gam_lam = expo_likelihood(1, gam_lam);
    p_gam_lam_prop = expo_likelihood(1, gam_lam_proposal);
    
    
    A1 = p_gam_lam_prop / p_gam_lam; 
    A2 = (gamma(gam_lam)/gamma(gam_lam_proposal))^p ;
    A3 = ((2*tau_prior_gamma_d^2)^(-p) * prod(1./diag(V_L)))^(gam_lam_proposal - gam_lam); 
    if (A3 == Inf); A3 = 1e+18;end;

    acceptance_probability_lambda_gamma_vec(irep) = min(1, A1*A2*A3 ); % A1*A2*A3 is a probability
    if binornd(1, acceptance_probability_lambda_gamma_vec(irep))  == 1
        if_accepted(irep) = 1;
        tau_prior_gamma_c = gam_lam_proposal;
    end
    
    e_star = 2 + p * tau_prior_gamma_c;
    f_star = M/(2 * tau_prior_gamma_c) + 1/2*(sum(1./diag(V_L)));
    tau_prior_gamma_d = sqrt(1/(Draw_iGamma(e_star, f_star))); % page 181

    % Now that we have the new estimate of tau2, update V_L (the prior
    % covariance matrix of beta)
    V_L = diag(tau2);
    
    % 4. Update sigma2 from Inverse Gamma
    c1 = (N-1+p)/2;
    PSI = (y-x*beta)'*(y-x*beta);
    c2 = 0.5*PSI + 0.5*(beta'/V_L)*beta;
    sigma2 = Draw_iGamma(c1,c2);
    
    % Save draws
    if irep > nburn
         beta_draws(irep-nburn,:) = beta;
         tau2_draws(irep-nburn,:) = tau2;
         sigma2_draws(irep-nburn,:) = sigma2;
    end
    beta_draws(irep,:) = beta;
    tau2_draws(irep,:) = tau2;
    sigma2_draws(irep,:) = sigma2;
  
end

model.beta = beta_draws;
model.sigma = sigma2_draws;

ols_beta = inv(x'*x)'*(x'*y);
mybeta = mean(beta_draws);
regression_accuracy(y', mybeta * x', 'mse')
regression_accuracy(y', ols_beta' * x', 'mse')
%regression_accuracy(y', beta_ls' * x', 'mse')

