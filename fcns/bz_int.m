function alpha_int = bz_int(alpha,x0,s_max)

if nargin == 2
    s_max = 1;
end

[n,m] = size(alpha);        % make sure alpha is a row vector
if n > m
    alpha = alpha';
end

M = length(alpha);
AA = zeros(M+1,M+1);

for ii = 1:M
    AA(ii,ii:ii+1) = [-1 1];
end

AA = M/s_max*AA;
AA(M+1,1) = 1;


alpha_int = (AA\[alpha x0]')';


