function Pri_up = Bay_2nd_mom(pho_SLD, models, L, n_modes, n_sam, n_max, n_iter, varargin)
% Importance sampling for getting the updated mean and covariance of the
% posterior.
%%%%%%%%%%%%%%%%% Parser %%%%%%%%%%%%%%%%%%%%%%%%%%
FWHM = 2.354820;
%seed = 12;
%rng(seed);


%defaultmin_sp = FWHM/2*0.02;
defaultmin_sp = 0;
defaultposCov = [];
defaultmax_r = 0.35;
defaultHG_flag = 0;
defaultcv_rate = [1,1,1];
defaultbri_known = 0;
defaultaperture = [0,0,1];
defaultU = 1;

p = inputParser;

addRequired(p,'pho_SLD');
addRequired(p,'models');
addRequired(p,'L');
addRequired(p,'n_modes');
addRequired(p,'n_sam');
addRequired(p,'n_max');


addOptional(p,'min_sp',defaultmin_sp);
addOptional(p,'posCov',defaultposCov);
addOptional(p,'max_r',defaultmax_r);
addOptional(p,'HG_flag',defaultHG_flag)
addOptional(p,'cv_rate',defaultcv_rate);
addOptional(p,'bri_known',defaultbri_known);
addOptional(p,'aperture',defaultaperture);
addOptional(p,'U',defaultU);


parse(p, pho_SLD, models, L, n_modes, n_sam, n_max, varargin{:});

pho_SLD = p.Results.pho_SLD;
models = p.Results.models;
L = p.Results.L;
n_modes = p.Results.n_modes;
n_sam = p.Results.n_sam;
n_max = p.Results.n_max;

min_sp = p.Results.min_sp;
posCov = p.Results.posCov;
max_r = p.Results.max_r;
HG_flag = p.Results.HG_flag;
cv_rate = p.Results.cv_rate;
bri_known = p.Results.bri_known;
aperture = p.Results.aperture;
U = p.Results.U;

perb_prob = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pri_up = zeros(size(models));

likelihood = zeros(size(models,3),1);

if n_iter == 1
    
    %n_sam = n_sam*2;
    
end

if ~HG_flag
    
    [V_L, D_L] = eig(L);
    
    D_L = real(D_L);
    
    [V_L, ~] = sortem(V_L, D_L);
    
else
    
    
    V_L = eye(size(L));
    
end

for i = 1:size(L,1)+1 % consider also the residual photons
    
    n(i) = sum(pho_SLD == i);
    
end

n_fac_max = 170;
% compute the log of the factorial prefactor in the multinomial log(N!/n1!...nM!)
while 1
    
    NN = log_factorial_app(sum(n),n_fac_max) - sum(log_factorial_app(n,n_fac_max));
    
    if ~isnan(NN) && NN ~= inf
        
        break;
        
    else
        
        n_fac_max = n_fac_max -1;
        
    end
    
end


for i = 1:size(models,3)
    
    err_sam_count = 0;
    
    mp = models(:,:,i) % pick a model
    
    mp = mp(1:nnz(mp(:,1)), :);
    
    pri_up_temp = zeros(size(mp));
    
    sam = zeros(size(mp,1), 3, n_sam);
    
    if size(mp, 1) == 1
        
        sam(:,1,:) = ones([1,1,n_sam]);
        
    else
        
        if bri_known
            
            sam(:,1,:) = ones([size(mp,1),1,n_sam])/size(mp,1);
            
        else
            
            sam(:,1,:) = permute( drchrnd((mp(:,2))', n_sam), [2,3,1] );
            
        end
        
        
        
    end
    
    
    % generation of samples
    
    if ~isempty(posCov)
        
        % Using covarience matrix
        
        if rand() <= perb_prob
            
            Cov_ratio = max(1, exp(log(cv_rate(2) - (cv_rate(2)- cv_rate(1))/cv_rate(3)*(n_iter-1) ) ) ) ;
            
            posCov{i} = posCov{i}*Cov_ratio;
            
        end
        
        if n_iter <=10 && max(eig(posCov{i}))< 1e-4
            
            posCov{i} = posCov{i}/max(eig(posCov{i}))*1e-3;
            
        end
        
        posCov_new{i} = 0;
        
        pos = mp(:,3:4);
        
        sam_temp = mvnrnd(pos(:), real(posCov{i}), n_sam );
        
        sam(:,2,:) = permute( sam_temp(:, 1:size(sam_temp,2)/2)  ,[2,3,1] );
        
        sam(:,3,:) = permute( sam_temp(:, size(sam_temp,2)/2+1:size(sam_temp,2) )  ,[2,3,1] );
        
    else
        
        % Using variance
        
        for j = 1:size(mp,1)
            
            sam(j,2,:) =  normrnd(mp(j,3), sqrt(mp(j,5)), [1,1,n_sam]);
            
            sam(j,3,:) =  normrnd(mp(j,4), sqrt(mp(j,6)), [1,1,n_sam]);
            
        end
        
        
    end
    
    % integration
    
    
    % calculation of log likelihood
    log_like = log_likelihood(n_modes, sam, V_L, n, aperture, U);
    log_like = real( log_like + NN );
    
    % Rescaling the likelihood
    
    cond1 = sum(exp(log_like - log(n_sam))) == inf;
    cond2 = any(isnan(exp(log_like - log(n_sam))));
    cond3 = sum(exp(log_like - log(n_sam))) == 0;
    
    if cond1 || cond2 || cond3
        
        [like_temp, like_scale(i)] = reSacle_like(log_like, n_sam, 1);
        
    else
        
        like_temp = exp(log_like - log(n_sam));
        like_scale(i) = 0;
        
    end
    
    sam(:,:,like_temp == 0) = [];
    like_temp(like_temp == 0) = [];
    
    % prevent out of bound
    out_bound = any( sqrt(sum(sam(:,2:3,:).^2,2)) > max_r);
    
    sam(:,:,out_bound) = [];
    like_temp(out_bound) = [];
    
    n_sam_eff = numel(like_temp);
    
    % update brightness
    
    if bri_known
        
        pri_up_temp(:,1) = ones(size(mp,1),1)/size(mp,1);
        
    else
        
        pri_up_temp(:,1) = sum( repmat( permute(like_temp,[1,3,2]), [size(sam,1), 1, 1]) ...
            .*sam(:,1,:)/sum(like_temp),3);
        
        pri_up_temp(:,1) = pri_up_temp(:,1)/sum(pri_up_temp(:,1));
        
    end
    
    
    
    % update position
    pri_up_temp(:,3:4) = sum( repmat( permute(like_temp,[1,3,2]), [size(sam,1), 2, 1]) ...
        .*sam(:,2:3,:)/sum(like_temp),3);
    
    
    pri_up_temp(:,2) = ( pri_up_temp(:,1) )*( sum(mp(:,2)) + 5*n_sam_eff/n_sam - size(mp,1)) + 1;
    
    % update SD
    
    if isempty(posCov)
        
        pri_up_temp(:,5:6) = sum( repmat( permute(like_temp,[1,3,2]), [size(sam,1), 2, 1]) ...
            .*( sam(:,2:3,:) -repmat(pri_up_temp(:,3:4),[1,1,n_sam_eff]) ).^2 ...
            /sum(like_temp),3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add randomness to the covariance to prevent position locking %%% NICO'S UPDATE%%%%
        pri_up_temp(:,5:6) =  pri_up_temp(:,5:6) + rand(size(pri_up_temp,1),2)*1e-6/n_iter;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        
        pos = reshape(pri_up_temp(:,3:4), [2*size(pri_up_temp,1), 1]);
        
        for pk = 1:n_sam_eff
            
            sam_temp = sam(:,2:3,pk);
            
            temp_sig = like_temp(pk)*( sam_temp(:) - pos(:) )*( sam_temp(:) - pos(:) )'/sum(like_temp);
            
            posCov_new{i} = posCov_new{i} + temp_sig;
            
        end
        
        if all(all(~isnan(posCov_new{i}))) && all(all(posCov_new{i}~=inf)) && all(all(posCov_new{i}~=-inf))
            
            [Cov_V, Cov_D] = eig(posCov_new{i});
            
            posCov_new{i} = Cov_V*( (real(Cov_D)>=0).*Cov_D + (real(Cov_D)<0).*zeros(size(Cov_D)) )*Cov_V';
            
        else
            
            posCov_new{i} = posCov{i};
            
        end
        
        pri_up_temp(:,5:6) = double(reshape(diag(posCov_new{i}), [size(mp,1), 2]));
        
    end
    
    pri_up(:,:,i) = zero_padding( pri_up_temp, n_max );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    post_sam = zeros(size(pri_up_temp,1), 3, n_sam);
    
    if size(pri_up_temp, 1) == 1
        
        post_sam(:,1,:) = ones([1,1,n_sam]);
        
    else
        
        if bri_known
            
            post_sam(:,1,:) = ones([size(mp,1),1,n_sam])/size(mp,1);
            
        else
            
            post_sam(:,1,:) = permute( drchrnd((pri_up_temp(:,2))', n_sam), [2,3,1] );
            
        end
        
    end
    
    for j = 1:size(pri_up_temp,1)
        
        post_sam(j,2,:) =  normrnd(pri_up_temp(j,3), sqrt(pri_up_temp(j,5)), [1,1,n_sam]);
        
        post_sam(j,3,:) =  normrnd(pri_up_temp(j,4), sqrt(pri_up_temp(j,6)), [1,1,n_sam]);
        
    end
    
    post_log_like = log_likelihood(n_modes, post_sam, V_L, n, aperture, U);
    
    post_log_like = post_log_like + NN;
    
    % prevent out of bound
    out_bound = any( sqrt(sum(post_sam(:,2:3,:).^2,2)) > max_r);
    
    post_sam(:,:,out_bound) = [];

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    likelihood(i) = sum(post_log_like)- log(n_sam);
    
    clear sam
    
end

Pri_up{1} = pri_up;
Pri_up{2} = likelihood;

if ~isempty(posCov)
    
    Pri_up{3} = posCov_new;
    
else
    
    Pri_up{3} = [];
    
end

end






