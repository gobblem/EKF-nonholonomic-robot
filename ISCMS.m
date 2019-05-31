function Sigma_Inv_Sparse = ISCMS( S , lambda , rho , plotting )
    [ N , M ] = size(S);
    if N~=M
        error('S is not square')
    end
    
    if sum(min(0,eig(S)))~=0
        error('S has negative eigenvalues')
    end
    if sum(sum(S - S.'))~=0
        error('S is not symmetrix')
    end

    
    Z = zeros(N); U = zeros(N);
    admm_iterations = 20;
    hist = zeros(admm_iterations,1);
    
    % start ADMM algorithm
    for q = 1:admm_iterations

        % X-update
        [ Q , V ] = eig(rho*(Z - U) - S);
        Xhat = diag(V);
        Xhat = (Xhat + sqrt(Xhat.^2+4*rho))/(2*rho);
        X = Q*diag(Xhat)*Q.';

        % Z-update
        Z = max(0,X+U-lambda/rho) - max(0,-X-U-lambda/rho);

        % U-update
        U = U+X-Z;
        
        hist(q) = norm(S*X);
        
    end
    
    if strcmp(plotting,'plot')
        plot(hist)
        title('norm of difference between identity and SX')
    end

    Sigma_Inv_Sparse = X;
end