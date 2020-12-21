function M = simple_elementwise(M)
    for i=1:size(M,1)
        for j=1:size(M,2)
            M(i,j) = simplify(M(i,j)) ;
        end
    end
end