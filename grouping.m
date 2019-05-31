function Groups = grouping( Sc , S )

    N = size(Sc);

    G = [];
    Groups = [];
    group = struct();
    group.members = [];
    q = 1;
    for j = N:-1:1
        % checks if the particular node has already been grouped
        if sum(find(G==j))
            continue
        end
        % if not already grouped, it creates a new group, and adds this node
        G = [ G , j ];
        Groups = [ Groups ; group ];
        for k = 1:j
            if Sc(k,j)~=0
                G = [ G , k ];
                Groups(q).members = [ Groups(q).members , k ];
            end
        end
        q = q + 1;
    end


end