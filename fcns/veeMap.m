
function out = veeMap(in)

out = zeros(3,1);

out(1) = -in(2,3);
out(2) = in(1,3);
out(3) = -in(1,2);

end