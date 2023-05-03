function HG_proj = HG_projection_GMM(n_modes, est, aperture, U)
% These are the expansion coeffs of the shifted PSF (to the source
% positions in est) into the tensor product HG modes for a gaussian mixture
% model of the aperture:
%   aperture - sub aperture positions and radii
%   U        - mixing unitary matrix 



n_ap = size(aperture,1);
[pj,qj,uj] = Indices_HG_GMMAperture(n_modes,n_ap); 
%est = est(1:nnz(est(:,1)),:);

% gpu setup
%{
if gpuDeviceCount("available") && size(est,3)>1
    est = gpuArray(est);
    aperture = gpuArray(aperture);
    U = gpuArray(U);
end
%} 


% parameters
s_b = reshape(est(:,1,:),[size(est,1),1,size(est,3)]);          % source brightnesses 
bar_mu = reshape(est(:,2:3,:),[size(est,1),2,1,size(est,3)]);   % source positions

% multi-aperture phase
Phi_xy = exp(1i*pagemtimes(aperture(:,1:2),'none', bar_mu,'transpose'));
%Phi_xy = exp(1i*aperture(:,1:2)*bar_mu');
B = 1/sqrt(n_ap) * pagemtimes(U , Phi_xy);


% single-aperture HG terms
pq = reshape([pj;qj],[1,2,numel(pj)]);            % HG powers for x and y
HG = prod(bar_mu.^pq .* exp(-bar_mu.^2/2)...
         ./ sqrt(factorial(pq)),2);



%HG_multiap = B(uj,:)'.*HG;
HG_multiap = permute(conj(B(uj,:,:,:)),[2,3,1,4]) .*HG;

HG_proj = gather(permute(HG_multiap,[1,3,4,2]));




%{
% ORIGINAL LINEAR INDEXING IMPLEMENTATION FOR MULTI_APERTURE - RIGHT BUT
SLOW
for count = 1:numel(uj)
    HG_proj(count).ind_u = uj(count);
    HG_proj(count).ind_x = pj(count);
    HG_proj(count).ind_y = qj(count);
    
    bar_mu  = est(:,2:3);              % mean of gaussian posteriors for each source
    s_b = est(:,1);                                    % brightness 
    
    
    HG_proj(count).xy = (bar_mu).^(repmat([HG_proj(count).ind_x, HG_proj(count).ind_y], [size(est,1),1])).*exp(-bar_mu.^2/2) ...
                               ./repmat(sqrt([factorial(HG_proj(count).ind_x),factorial(HG_proj(count).ind_y)]), [size(est,1),1]);
    
    % bar_mu .^ [p,q] .* exp(-bar_mu.^2 / 2) ./ sqrt([p!,q!])  
    % A = nx2
    % R = mx2
    % A*R' = nxm (each column is a unique source)
    % U = nxn
    % U * exp(A*R') = nxm
    Phi_xy = exp(1i*aperture(:,1:2)*bar_mu');
    
    B = 1/sqrt(n_ap) * U(HG_proj(count).ind_u,:) * Phi_xy;
    
    HG_proj(count).proj = conj(B).' .* prod(HG_proj(count).xy, 2);
    
    HG_proj(count).total = sum( s_b .* HG_proj(count).proj, 1); 
    
    %HG_proj(count).proj = conj(phase_fn(bar_mu',aperture,U)).' .* prod(HG_proj(count).xy, 2);
    
    %HG_proj(count).total = sum(s_b .* prod( conj(phase_fn(bar_mu', aperture,U)).' .* HG_proj(count).xy, 2).^2, 1);    % multiply in source brighnesses
    
    
    %HG_proj(count).xy = (est(:,2:3)).^(repmat([HG_proj(count).ind_x, HG_proj(count).ind_y], [size(est,1),1])).*exp(-est(:,2:3).^2/2) ...
    %                           ./repmat(sqrt([factorial(HG_proj(count).ind_x),factorial(HG_proj(count).ind_y)]), [size(est,1),1]);

    %HG_proj(count).proj = prod(conj(phase_fn(est,aperture,U)).' .* HG_proj(count).xy, 2);

    %HG_proj(count).total = sum(est(:,1) .* prod( conj(phase_fn(est,aperture,U)).' .* HG_proj(count).xy, 2).^2, 1);    % multiply in source brighnesses    
end
%}



end