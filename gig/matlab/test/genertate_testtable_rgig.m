function [MC_res, numerical_res] = genertate_testtable_rgig(p, a, b, sim)
%%
%
%   Test the simulator 
%   First compare 
%   Ex, Ex^{-1}, Elog(x) vs numerical values
%
%   in:
%    p   - (nx1) param
%    a   - (nx1) param
%    b   - (nx1) param
%    sim - int  number MC simulations
%  
%   out:
%   MC_res
%   (1) - (nx1)  MC Ex
%   (2) - (nx1)  MC Ex^-1
%   (3) - (nx1)  MC Elog(x)
%
%   numerical_res
%   (1) - (nx1)  num Ex
%   (2) - (nx1)  num Ex^-1
%   (3) - (nx1)  num Elog(x)
%%
n = length(p);
MC_res = zeros(n,3);
numerical_res = zeros(n,3);
[numerical_res(:,2), numerical_res(:,1) ,numerical_res(:,3)] = Egig(p , a, b);
for i = 1:sim
    gam_sim = rgig_par(p, a, b,ceil(10^6*rand(6,1)) );
    
    MC_res(:,1)  = MC_res(:,1)  + gam_sim;
    MC_res(:,2)  = MC_res(:,2)  + 1./gam_sim;
    MC_res(:,3)  = MC_res(:,3)  + log(gam_sim); 
end
MC_res = MC_res/sim;
%computing probablilites numerically

% for i = 1:n
%     % [0,Ex/2]
%     quadl(@(x)gig_den(x,p(i),a(i),b(i)) ,0, numerical_res(i,1));
%     % [0,Ex]
%     quadl(@(x)gig_den(x,p(i),a(i),b(i)) ,0, numerical_res(i,1));
%     % [0,2*Ex]
%     quadl(@(x)gig_den(x,p(i),a(i),b(i)) ,0, numerical_res(i,1));  
% end
%