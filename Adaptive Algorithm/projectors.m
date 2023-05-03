function projectors = projectors(tp_temp)
    % This function computes the eigenvectors of the Bayesian Quantum 
    % Cramer-Rao Lower Bound matrix (Sigma_Q). The minimum eigenvector of
    % this matrix ends up being the weights to the linear combination of
    % estimation parameters that comprise the joint parameter.
    n_now = size(tp_temp,1);
    
    L(:,:,:) = SLD_3N(tp_temp(:,1),tp_temp(:,3:4),20); % tgammas is the SLD operators
     
     for i = 1:3*n_now
         for j = 1:3*n_now
             
             H(i,j) = trace( L(:,:,3*n_now+1)*( L(:,:,i)*L(:,:,j) ...
                           + L(:,:,j)*L(:,:,i) )/2 ); % tgammas is the QFI matrix elements
             
         end
     end
     
     %H(:,:,k) = real(H(:,:,k));
     
     [temp_V,temp_D] = eig(H); % the eigenvalues and eigenvectors of QFI matries
     [temp_V,temp_D] = sortem(temp_V,temp_D);
     
     projectors = temp_V(:,1);


end