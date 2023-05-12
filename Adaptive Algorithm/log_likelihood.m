function log_like = log_likelihood(n_modes, sam, V_L, n, aperture, U)
  
    % gpu setup
    %{
    if gpuDeviceCount("available")
        gpuDevice(1)
        sam = gpuArray(sam);
        V_L = gpuArray(V_L);
    end
    %}
    
    if gpuDeviceCount("available")
        max_batch_size = floor((3e9/8) / (numel(V_L) * prod(size(sam,1,2))));
    else
        max_batch_size = floor((7e9/8) / (numel(V_L) * prod(size(sam,1,2))));
    end

    n_batches = floor(size(sam,3)/max_batch_size);
    rem_batch_size = rem(size(sam,3),max_batch_size);
    
    k = 1;
    for b = 1:n_batches+1
        
        % check if we are at the remainder batch and change the batch size
        if b == n_batches + 1
            max_batch_size = rem_batch_size;
        end
        
        % turn the samples of the batch into a gpu array if possible
        if gpuDeviceCount("available")
            sam_batch = gpuArray(sam(:,:,k:(k+max_batch_size-1)));
        else
            sam_batch = sam(:,:,k:(k+max_batch_size-1));
        end 
        

        % compute density operator for each sample
        rho_sam_batch =  rho_HG_GMM(n_modes,sam_batch,aperture,U);
        
        
        
        % get modal probabilities for each sample given the measurement matrix V_L
        prob_sam_batch = real(sum(conj(V_L).* pagemtimes(rho_sam_batch,V_L),1));
        prob_conj_batch = 1 - sum(prob_sam_batch,2);
        prob_sam(1,:,k:(k+max_batch_size-1)) = cat(2,prob_conj_batch.*double(prob_conj_batch>0),prob_sam_batch);
        
        % remove negative probabilities
        prob_sam(prob_sam < 0) = 0;
        
        % increment batch index
        k = k + max_batch_size;
    end
    
    %{
    % compute density operator for each sample
    %rho_sam = rho_HG(n_modes, sam);
    rho_sam = rho_HG_GMM(n_modes,sam,aperture,U);
    
    % get modal probabilities for each sample given the measurement matrix
    % V_L
    prob_sam = real(sum(conj(V_L).* pagemtimes(rho_sam,V_L),1));
    prob_conj = 1 - sum(prob_sam,2);    
    prob_sam = cat(2,prob_conj.*double(prob_conj>0),prob_sam);
    %}
    
    % compute log probabilties given the measuremenet vector n
    log_prob = n.* log(prob_sam + (prob_sam == 0));
    
    % compute the log likelihood of each sample
    log_like = gather(squeeze(sum(log_prob))');
    
            
    
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