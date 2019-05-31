function p = distributed_KF( p , rho  )
%A , B , C , D , u_t , mu_tm_tm , sigma_tm_tm , y , Q , R )

    % Breaks out the relevant arguments from the input data structure
    A = p.A;
    B = p.B;
    D = p.D;
    u = p.u;
    
    K = p.K;
    
    mu_tm_tm = p.mu_prior;
    sigma_tm_tm = p.sigma_prior;
    R = p.R_global;
    
    % Initialization
    N = length(mu_tm_tm);
    T1 = p.L;
    T2 = 1;

    p.mu_hist = zeros(N,K);
    p.sigma_hist = zeros(N,N,K);
    
    % Runs KF over entire history of measurements
    for t = 1:p.K
        X_mean = zeros(N,1);
        X = zeros(N,T1+T2);
        Z = zeros(N,T1+T2);

        % ADMM loop to solve for mean: 
        for it = 1:50
            mean_index = zeros(N,1);
            % This loop can be parallelized
            for k = 1:T1
                Q = p.partition_struct(k).Q;
                y = p.partition_struct(k).y(:,t);
                L = length(y);
                mask = p.partition_struct(k).mask;
                C = p.partition_struct(k).C(:,mask);
                mean_index = mean_index + double(mask).';
                X(mask,k) = ( C.'*(Q\C)+eye(L)*rho )\( C.'*(Q\y)+rho*X_mean(mask) - Z(mask,k) );
            end
            % This loop can be parallelized
            for k = 1:T2
                sigmahat = A*sigma_tm_tm*A.' + R;
                muhat = A*mu_tm_tm+B*u(:,t);
                X(:,T1+k) = (  sigmahat + eye(N)*rho )\( sigmahat*muhat - Z(:,T1+k) + rho*X_mean );
                mean_index = mean_index + T2;
            end

            % This mean calculation can be parallelized
            X_mean = sum(X./repmat(mean_index,1,T1+T2),2);
            
            for k = 1:T1
                mask = p.partition_struct(k).mask;
                Z(mask,k) = Z(mask,k) + rho*(X(mask,k)-X_mean(mask)); 
            end
            for k = 1:T2
                Z(:,T1+k) = Z(:,T1+k) + rho*(X(:,T1+k) - X_mean);
            end


        end

        % Loops to solve for new covariance matrix: This can be done in
        % parallel with finding the mean
        sigma_t_tm = A*sigma_tm_tm*A.' + R;
        sigma_t_t = ( eye(N) - sigma_t_tm*p.C.'*((p.C*sigma_t_tm*p.C.' + p.Q_global)\p.C))*sigma_t_tm;
        
        p.sigma_hist(:,:,t) = sigma_t_t;
        
        p.mu_hist(:,t) = X_mean;
        
        mu_tm_tm = X_mean;
        
        sigma_tm_tm = sigma_t_t;
        

    end
    
    p.mu = p.mu_hist;
    p.sigma = p.sigma_hist;
    
end