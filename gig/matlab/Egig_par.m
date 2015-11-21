function [EG,EiG, ElogG] = Egig_par(p, a, b)
%%WARNING 
% 
%	Calculates the expcation from GIG distribution has the pdf
%   f(x|p, a, b) = \frac{(a/b)^{p/2}}{2K_{p}(\sqrt{ab})} x^{p-1} e^{-(ax + bx^{-1})/2}
%	where K_p is modifed Bessel function of second (third?) kind.
%	the parameters should be in vector format
% 
%   p	 -(nx1) parameter in f(x)
%	a	 -(nx1) parameter in f(x)
%	b	 -(nx1) parameter in f(x)
%
%
% Returns:
%
%
%	Eg	- nx1 E[G]
%	EiG	- nx1 E[G^-1]
%	ElogG	- nx1 E[log(G)]
%
