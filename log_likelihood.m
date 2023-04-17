function log_likelihood = log_likelihood(n_modes, sam, V_L, n)

    n_sam = size(sam, 3);

    for k = 1:n_sam
       
        %k
        
        rho_temp = rho_HG(n_modes, sam(:,:,k));
        
        prob_temp = real( diag(V_L'*rho_temp*V_L)' );
        
        prob_temp = [max(0, 1 - sum(prob_temp)), prob_temp]; % consider also the residual photons
        
        log_prob_temp = n.*log(prob_temp + (prob_temp == 0));
        
        log_likelihood(k) = sum(log_prob_temp);
        
    end


end