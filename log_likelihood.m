function log_likelihood = log_likelihood(n_modes, sam, V_L, n)
  
    
    % compute density operator for each sample
    rho_sam = rho_HG(n_modes, sam);
    
    % get modal probabilities for each sample given the measurement matrix
    % V_L
    prob_sam = real(sum(conj(V_L).* pagemtimes(rho_sam,V_L),1));
    prob_conj = 1 - sum(prob_sam,2);
    prob_sam(:,end+1,:) = prob_conj.*double(prob_conj>0); 
    
    
    % compute log probabilties given the measuremenet vector n
    log_prob = n.* log(prob_sam + (prob_sam == 0));
    
    % compute the log likelihood of each sample
    log_likelihood = squeeze(sum(log_prob))';
    
    
    
   %{
    %OLD IMPLEMENTATION
    n_sam = size(sam, 3);

    for k = 1:n_sam
       
        %k
        
        rho_temp = rho_HG(n_modes, sam(:,:,k));
        
        prob_temp = real( diag(V_L'*rho_temp*V_L)' );
        
        prob_temp = [max(0, 1 - sum(prob_temp)), prob_temp]; % consider also the residual photons
        
        log_prob_temp = n.*log(prob_temp + (prob_temp == 0));
        
        log_likelihood(k) = sum(log_prob_temp);
        
    end
   %}

end