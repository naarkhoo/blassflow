function [Eig,Eg, Elogg] = Egig(p, a, b)
%%
%
%   Egig(p, a, b)
%
%   computing expctation of generalised inverse gaussian distribution
%
%%

sqrt_ab = sqrt(a.*b);

Kp = besselk(p, sqrt_ab);
Kp1 = besselk(p + 1, sqrt_ab);

Eg = (b ./ a).^(1/2) .* Kp1 ./ Kp;

Kp1 = besselk(p - 1, sqrt_ab);
Eig = (b ./ a).^(-1/2) .* Kp1 ./ Kp;
e =1e-6;

Kpe = besselk(p + e, sqrt_ab);
Elogg = log(b./a)/2 + (1/e) * (Kpe ./ Kp - 1);