function EA = fcn_X2EA(X)

vR = X(7:15);
R = reshape(vR,[3,3]);
EA = veeMap(logm(R))';


end
