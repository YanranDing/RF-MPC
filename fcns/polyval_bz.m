% Function to evaluate bezier polynomials
% Inputs: Alpha - Bezeir coefficients (alpha_0 ... alpha_M)
%         s - s parameter. Range [0 1]
% Outputs: b = sum(k=0 to m)[ alpha_k * M!/(k!(M-k)!) s^k (1-s)^(M-k)]
%factorial(M)/(factorial(k)*factorial(M-k))
function b = polyval_bz(alpha, s)
    b = zeros(size(s)) ;
    M = size(alpha,2)-1 ;  % length(alpha) = M+1
    switch M
        case 2
            c = [1 2 1];
        case 3
             c = [1 3 3 1];
        case 4
            c = [1 4 6 4 1];
        case 5
            c = [1 5 10 10 5 1];
        case 6
            c = [1 6 15 20 15 6 1];
        case 7
            c = [1 7 21 35 35 21 7 1];
        case 8
            c = [1 8 28 56 70 56 28 8 1];
        case 9
            c = [1 9 36 84 126 126 84 36 9 1];
        case 10
            c = [1 10 45 120 210 252 210 120 45 10 1];
        otherwise
            c = 0;
    end
        
    
        
    for k = 0:M
        
        b = b + alpha(:,k+1) .* factorial(M)/(factorial(k)*factorial(M-k)) .* s.^k .* (1-s).^(M-k) ;
    end
    
%    b = alpha(:,1) .* c(1) .* (1-s).^(M)+alpha(:,2) .* c(2) .* s.^1 .* (1-s).^(M-1) + alpha(:,2+1) .* c(2+1) .* s.^2 .* (1-s).^(M-2) + alpha(:,3+1) .* c(3+1) .* s.^3 .* (1-s).^(M-3) + alpha(:,4+1) .* c(4+1) .* s.^4 .* (1-s).^(M-4) + alpha(:,5+1) .* c(5+1) .* s.^5 .* (1-s).^(M-5);
end