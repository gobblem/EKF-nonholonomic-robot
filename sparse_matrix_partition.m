function fstruct = sparse_matrix_partition( G , S )

    fstruct = struct();

    L = length(G);
    for k = 1:L
        n = length(G(k).members);
        fstruct(k).Q = zeros(n);
        m = sort(G(k).members);
        fstruct(k).m = m;
        for i = 1:n
            for j = 1:n
                fstruct(k).Q(i,j) = S(m(i),m(j));
            end
        end
    end

end