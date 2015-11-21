function x = rgig_par(p, a, b, seed)
%%WARNING 
% 
%	Samples from GIG distribution has the pdf
%   f(x|p, a, b) = \frac{(a/b)^{p/2}}{2K_{p}(\sqrt{ab})} x^{p-1} e^{-(ax + bx^{-1})/2}
%	where K_p is modifed Bessel function of second (third?) kind.
%	the parameters should be in vector format
% 
%   Algorithm created from:
%   Generating Generalized Inverse Gaussian Random Variates, H??rmann and Leydold, 2013
%
%   p	 -(nx1) parameter in f(x)
%	a	 -(nx1) parameter in f(x)
%	b	 -(nx1) parameter in f(x)
%	seed -(6x1) int (the seed)
%
%
x = gig_matlab(p, a, b);