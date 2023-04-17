function HG_proj = HG_projection_GMM(n_mode, est)
% These are the expansion coeffs of the shifted PSF (to the source
% positions in est) into the tensor product HG modes for a gaussian mixture
% model of the aperture.


% aperture.mat must contain
%   1) sub aperture positions and radii
%   2) mixing unitary matrix 
load 'aperture.mat';

est = est(1:nnz(est(:,1)),:);

%{
count = 1;

for k = 1:n_apertures
    for i = 1:n_mode
        for j = 1:i
            HG_proj(count).ind_k = k;
            HG_proj(count).ind_x = i-j;
            HG_proj(count).ind_y = i-HG_proj(count).ind_x-1;

            HG_proj(count).xy = (est(:,2:3)).^(repmat([HG_proj(count).ind_x, HG_proj(count).ind_y], [size(est,1),1])).*exp(-est(:,2:3).^2/2) ...
                               ./repmat(sqrt([factorial(HG_proj(count).ind_x),factorial(HG_proj(count).ind_y)]), [size(est,1),1]);

            HG_proj(count).proj = prod(conj(phase_fn(est,aperture,U)).* HG_proj(count).xy, 2);

            HG_proj(count).total = sum(est(:,1).*prod( conj(phase_fn(est,aperture,U)) .* HG_proj(count).xy, 2).^2, 1);    % multiply in source brighnesses

            count = count + 1;

        end
    end
end
%}

for count = 1:N_modes
    HG_proj(count).ind_u = uj(count);
    HG_proj(count).ind_x = pj(count);
    HG_proj(count).ind_y = qj(count);
    
    
    HG_proj(count).xy = (est(:,2:3)).^(repmat([HG_proj(count).ind_x, HG_proj(count).ind_y], [size(est,1),1])).*exp(-est(:,2:3).^2/2) ...
                               ./repmat(sqrt([factorial(HG_proj(count).ind_x),factorial(HG_proj(count).ind_y)]), [size(est,1),1]);

    HG_proj(count).proj = prod(conj(phase_fn(est,aperture,U)).' .* HG_proj(count).xy, 2);

    HG_proj(count).total = sum(est(:,1) .* prod( conj(phase_fn(est,aperture,U)).' .* HG_proj(count).xy, 2).^2, 1);    % multiply in source brighnesses    
end


end