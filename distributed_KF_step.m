function [ mu , sigma ] = distributed_KF_step( A , B , C , D , u , y , Q , R , K , mu_tm_tm , sigma_tm_tm )

    % Builds the data structure to pass to the solver
    p = struct();
    p.A = A;
    p.B = B;
    p.C = C;
    p.D = D;
    p.Q_global = Q;
    p.R_global = R;
    p.u = u;
    p.y = y;
    p.mu_prior = mu_tm_tm;
    p.sigma_prior = sigma_tm_tm;
    
    % This is the number of timesteps to run the KF for. This should match
    % the lengths of y and u
    p.K=K;

    % Finds sparse approximation to Q
    % parameters are ( Q , lambda , rho , 'plot' if you want an ADMM
    % converence plot or anything else otherwise. Tune lambda to adjust the
    % sparsity regularization parameter and rho to adjust the L2
    % regularization parameter to make the entries small.
    %---------Tune This------------
    S_sparse = ISCMS( Q , 2 , 1 , 'none');
    S_sparse = Q.*(abs(S_sparse)>1e-3);
    
    % Decomposes Q into disconnected parts
    Sc = chol(S_sparse);
    G = grouping(Sc);
    p.partition_struct = sparse_matrix_partition( G , Q );

    p.L = length(p.partition_struct);
    L = p.L;
    
    % Decomposes y and C and builds a mask for each decomposed part
    for k = 1:L
        p.partition_struct(k).y = p.y(p.partition_struct(k).m,:);
        p.partition_struct(k).C = p.C(p.partition_struct(k).m,:);
        p.partition_struct(k).mask = sum(p.partition_struct(k).C,1)~=0;
    end

    % Solves distributed KF
    % The parameter is rho which is the ADMM regularization parameter.
    % Changing this will affect the convergence rate of ADMM, so tune on a
    % test case with similar noise to a validation case
    %------Tune this------
    p = distributed_KF( p , 10 );

    mu = p.mu;
    sigma = p.sigma;

end