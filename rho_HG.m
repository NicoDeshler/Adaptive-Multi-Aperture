function rho_HG = rho_HG(n_modes, scene)

%    The 2 assignements below are only necessary if the number of sources is unknown apriori
%    scene = sortrows(scene, 'descend');
%    scene = scene(1:nnz(scene(:,1)), :);
    
%   HG_temp = HG_projection(n_modes, scene);
    HG_modes = HG_projection_GMM(n_modes,scene);
    %HG_temp = HG_projection_GMM(n_modes,scene);

    %HG_modes = [HG_temp.proj];
    
    % Density operator in HG mode basis
    %{
    rho = zeros(size(HG_modes,2));
    
    
    for j = 1:size(scene,1)
       
        rho = rho + scene(j,1)*(HG_modes(j,:))'*HG_modes(j,:); 
        
    end
    
    %}  
    
    rho_HG = pagemtimes( (scene(:,1,:).* HG_modes),'ctranspose',HG_modes,'none');
    
    % tr = arrayfun(@(k) real(trace(rho_HG(:,:,k))),1:size(rho_HG,3));
    %rho_HG = HG_modes'*diag(scene(:,1))*HG_modes;
end