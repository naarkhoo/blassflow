function f = gig_den(x, p, a , b)
sqrt_ab = sqrt(a.*b);

K_nom = (a./b).^(p/2);
K_den = 2* besselk(p, sqrt_ab);

f = (K_nom./K_den) .* x.^(p-1).*exp( - (a .* x + b./x)/2); 