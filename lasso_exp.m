% baysian lasso
% based on Korobilis, D  Hierarchical Shrinkage Priors for Dynamic Regressions with Many Predictors. International Journal of Forecasting 2013 implementation
% https://sites.google.com/site/dimitriskorobilis/matlab/bayesian-lasso

addpath('generators')
n_simulation = 10;

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
nsave = 2000;
nburn = 2000;
ntot = nsave + nburn;

beta_draws = zeros(nsave,p);
tau2_draws = zeros(nsave,p);
lambda2_draws = zeros(nsave,1);
sigma2_draws = zeros(nsave,1);
% ----------------Set priors
% lambda2 ~ Gamma(r,d)
r = 1;
delta = 1.78;

% Get OLS quanities from the full model (only if degrees of freedom allow this)
if N>p
    beta_OLS = inv(x'*x)'*(x'*y);
    SSE_OLS = (y - x*beta_OLS)'*(y - x*beta_OLS);
    sigma2_OLS = SSE_OLS/(T-(p-1));
end

% Initialize parameters
beta = 0*ones(p,1); %beta_OLS;
tau2 = 4*ones(p,1);
V_L = diag(tau2);
lambda2 = 0.1; %(p*sqrt(sigma2_OLS)/(sum(abs(beta_OLS)))).^2;
sigma2 = 0.1; %sigma2_OLS;

%==========================================================================
%====================| GIBBS ITERATIONS START HERE |=======================
disp('Now you are running the Bayesian Lasso with normal-exp prior')
for irep = 1:ntot
    if mod(irep,5000)==0
        disp(['Iteration: ' num2str(irep) ' of ' num2str(ntot) '. Please wait a few seconds...'])
    end
    
    beta_old = beta;
    
    % 1. Update beta from Normal
    A = inv(x'*x + inv(V_L));
    post_mean_beta = A*x'*y; %#ok<*MINV>
    post_var_beta = sigma2*A;
    beta = Draw_Normal(post_mean_beta,post_var_beta);
    
    % 2. Update tau2_j from Inverse Gaussian
    for j = 1:p
        a1 = (lambda2*sigma2)./(beta(j,1).^2);
        a2 = lambda2;
        tau2_inverse = Draw_IG(sqrt(a1),a2);
        tau2(j,1) = 1/tau2_inverse + 1e-15;
    end
    
    % Now that we have the new estimate of tau2, update V_L (the prior
    % covariance matrix of beta)
    V_L = diag(tau2);
    
    % 3. Update lambda2 from Gamma
    b1 = p + r;
    b2 = 0.5*sum(tau2) + delta;
    lambda2 = Draw_Gamma(b1,b2);
    
    % 4. Update sigma2 from Inverse Gamma
    c1 = (N-1+p)/2;
    PSI = (y-x*beta)'*(y-x*beta);
    c2 = 0.5*PSI + 0.5*(beta'/V_L)*beta;
    sigma2 = Draw_iGamma(c1,c2);
    
    % Save draws
    if irep > nburn
        beta_draws(irep-nburn,:) = beta;
        tau2_draws(irep-nburn,:) = tau2;
        lambda2_draws(irep-nburn,:) = lambda2;
        sigma2_draws(irep-nburn,:) = sigma2;
    end
    
end

ols_beta = inv(x'*x)'*(x'*y);
mybeta = mean(beta_draws);
regression_accuracy(y', mybeta * x', 'mse')
regression_accuracy(y', ols_beta' * x', 'mse')

