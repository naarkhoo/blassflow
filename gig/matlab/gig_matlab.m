function x = gig_matlab(p, a ,b)

n = length(a);
alpha = sqrt(a./b);
beta  = sqrt(a.*b);
p_abs = abs(p);
x = zeros(n, 1);
index1 = beta >1 | p_abs > 1;
if sum(index1) > 0
    x(index1) = algortihm_region1(p_abs(index1), beta(index1));
end
index2 = logical((beta <= min(1/.2, 2./3. * sqrt(1-p_abs))).*(1-index1));
if sum(index2)>0
    x(index2) = algorithm_region2(p_abs(index2),beta(index2));
end
index3 = logical(1 - index1 - index2);
if sum(index3) >0
    x(index3) = algortihm_region3(p_abs(index3),beta(index3));
end
index = p < 0;
x(index) = 1 ./ x(index);

x = x./alpha;


end

function x = algorithm_region2(p, beta)

two_d_beta = 2 ./ beta;
mode = beta;
mode = mode ./ (sqrt((1 - p).^2 + beta.^2 ) + (1 - p));
x0 = beta ./ (1 -p);
xs = max(x0, two_d_beta);

c1 = gig_propto(mode, p, two_d_beta);

A1 = c1 .* x0;

n = length(beta);

c2 = zeros(n, 1);
A2 = c2;

index_beta =  x0 < two_d_beta;
index_p    =  logical((p == 0).*index_beta);

c2(index_beta) = exp(-beta(index_beta));
x0_pow_p = x0(index_beta).^p(index_beta);
A2(index_beta) = c2(index_beta).* (two_d_beta(index_beta).^p(index_beta) - ...
                                   x0_pow_p) ./ p(index_beta);
A2(index_p) = c2(index_p) .* log(two_d_beta(index_p) ./ beta(index_p));

c3 = xs.^(p-1);
A3 = 2 * c3 .* exp( - xs ./ two_d_beta) ./beta;
A  = A1 + A2 + A3;
k = length(beta);
x = zeros(k,1);
index = 1:k;

p_0 = sum(p==0);

while( k> 0)
    U = rand(k,1);
    V = rand(k,1).*A;
    x_ = zeros(k,1);
    c4 = x_;
    index1 = V <= A1;

    x_(index1) = x0(index1) .* V(index1) ./ A1(index1);
    c4(index1) = c1(index1);

    index2 = (V <= A1 + A2).*(1-index1);
    index2 = logical(index2);
    V(index2) = V(index2) - A1(index2);

    index2_p = logical(index2 .* (p>0));
    if sum(index2_p) >0
        x_(index2_p) = (x0_pow_p(index2_p) + V(index2_p) .* p(index2_p) ./...
                       c2(index2_p)).^(1./p(index2_p));
    end
    if p_0> 0
        index2_np = logical(index2 .* (p==0));
        if sum(index2_np) > 0
            x_(index2_np) = beta(index2_np) .* exp(V(index2_np) .* exp(beta(index2_np)));

        end
    end
    c4(index2_p) = c2(index2_p) .* x_(index2_p).^(p(index2_p) - 1);
    index3 = logical(1 - index1 - index2);
    V3 = V(index3) - A1(index3) - A2(index3);
    twd = two_d_beta(index3);
    x_(index3) = - twd .* log( exp( - xs(index3) ./ twd) - V3 ./ (c3(index3).* twd));
    c4(index3) = c3(index3) .* exp( - x_(index3) ./ twd);
    
    index_acc = U .* c4 <= gig_propto(x_, p, two_d_beta);
    x(index(index_acc))   = x_(index_acc);
    x0(index_acc)         = [];
    A(index_acc)          = [];
    A1(index_acc)         = [];
    A2(index_acc)         = [];
    c2(index_acc)         = [];
    c3(index_acc)         = [];
    c1(index_acc)         = [];
    beta(index_acc)       = [];
    p(index_acc)          = [];
    xs(index_acc)         = [];
    x0_pow_p(index_acc)   = [];
    two_d_beta(index_acc) = [];
    index(index_acc)      = [];
    k = length(index);
end
end

function x = algortihm_region3(p, beta)

    two_d_beta = 2 ./ beta;
    mode = beta;
    mode = mode ./ (sqrt((1 - p).^2 + beta.^2 ) + (1 - p));
    
    x_p = 1 + p + sqrt( (1 + p).^2 + beta.^2);
    x_p = x_p./ beta;
    v_p = sqrt_gig_propto(mode, p, two_d_beta);
    u_p = x_p .* sqrt_gig_propto(x_p, p, two_d_beta);
    
    k = length(beta);
    x = zeros(k,1);
    index = 1:k;
    while( k> 0)
       U = rand(k,1) .* u_p;
       V = rand(k,1) .* v_p;
       x_ = U ./ V;
       index_acc = V <= sqrt_gig_propto(x_, p, two_d_beta);
       x(index(index_acc)) = x_(index_acc);
       mode(index_acc)       = [];
       p(index_acc)          = [];
       u_p(index_acc)        = [];
       v_p(index_acc)        = [];
       two_d_beta(index_acc) = [];
       index(index_acc)      = [];
       k = length(index); 
    end

end

function x = algortihm_region1(p, beta)
    two_d_beta = 2 ./ beta;
    mode = sqrt((p - 1).^2 + beta.^2 ) + (p - 1);
	mode = mode./beta;
    
    c1 = - two_d_beta .* (p + 1)  - mode;
	c2 =   two_d_beta .* (p - 1)  .* mode - 1;
    c3 = c2 - (c1.^2) / 3;
	c4 = 2 * (c1.^3) / 27 - c1 .* c2 / 3 + mode;
	c5 = acos(- c4./2 .* sqrt(- 27./ c3.^3 ));
 	c6 = sqrt(-4 * c3 / 3 );
	x_m = c6 .* cos( c5 / 3 + 4 * pi / 3   ) - c1 / 3;
	x_p = c6 .* cos( c5 / 3                ) - c1 / 3;
    
    % triangle to simulate
	u_m_div_v = (x_m - mode) .* sqrt_gig_ratio(x_m, mode, p, two_d_beta);
	u_p_div_v = (x_p - mode) .* sqrt_gig_ratio(x_p, mode, p, two_d_beta);
	u_p_m = u_p_div_v - u_m_div_v;
    
    
    k = length(beta);
    x = zeros(k,1);
    index = 1:k;
    while k > 0
    
        U = rand(k,1);
        V = rand(k,1);
        x_ = u_p_m .* (U ./ V) + (u_m_div_v ./ V) + mode;
        
        index_acc = V <= sqrt_gig_ratio(x_, mode, p, two_d_beta);
        index_acc(x_<0) = 0;
        x(index(index_acc)) = x_(index_acc);
        mode(index_acc)       = [];
        p(index_acc)          = [];
        two_d_beta(index_acc) = [];
        u_m_div_v(index_acc)  = [];
        u_p_m(index_acc)      = [];
        index(index_acc)      = [];
        k = length(index);
    end
    
end

function ret = sqrt_gig_ratio(x_in, m_in, p, two_d_beta)
    ret = (x_in ./ m_in).^((p - 1)/2);
	ret = ret.*exp( ( m_in + 1./m_in - x_in - 1./x_in  )./(2 * two_d_beta ) );
end

function ret = sqrt_gig_propto(x_in, p, two_d_beta)
   ret = exp( ((p-1)/2) .* log(x_in) - (x_in + 1 ./x_in) ./ (2 * two_d_beta));
end

function ret = gig_propto(x_in, p, two_d_beta)
   ret = exp( ((p-1)) .* log(x_in) - (x_in + 1 ./x_in) ./ ( two_d_beta));
end