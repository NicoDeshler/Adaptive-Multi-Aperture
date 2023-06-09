function AdaptiveArray_ErrorTest(array_id,num_workers)

    % add directories with all algorithm functions
    addpath('Adaptive Algorithm/')

    % get the job array id
    if ischar(array_id)
        array_id = str2double(array_id);
    end
       
    % make the data storing object
    DS = DSAdaptive_ErrorTest();

    % make the save directory
    mkdir(DS.save_dir)
    
    % get configuration indices
    [a,n,m,p,e1,e2] = ind2sub(DS.cfg_size,array_id);
    cfg_id = {a,n,m,p,e1,e2};
    
    
    % parfor configuration variables
    trials = DS.trials;
    aperture = DS.apertures{a};
    num_src = DS.num_src(n);
    min_sep_frac = DS.min_sep_frac(m);
    num_pho = DS.num_pho(p);
    dark_pho_lambda = DS.dark_pho_lambda(e1);
    phase_err_sigma = DS.phase_err_sigma(e2);
    max_order = DS.max_order;
    U = dftmtx(size(aperture,1))/sqrt(size(aperture,1));

    %------------ Rescale dimensions ---------------- %
    
    % Gaussian widths
    tilde_sigma = min(aperture(:,3))/3;  % stdev of gaussian sub-aperture
    tilde_sigma_eff = DS.R_max/3;        % stdev of the effective gaussian aperture
    sigma = 1/2/tilde_sigma;             % stdev of sub-aperture gaussian PSF
    sigma_eff = 1/2/tilde_sigma_eff;     % effective stdev of effective gaussian PSF
    
    % make units of the aperture coordinates equal to [(2 sigma)^-1] to
    % correspond to the fact that Kwan defines the image-space coordinates
    % to be [x/(2 sigma)].
    aperture = aperture ./ tilde_sigma;   
    % make the min_sep_frac of the effective rayleigh length scaled in
    % terms of the minimum separation fraction in terms of the sub-aperture
    % rayleigh length
    
    
    min_sep_frac = min_sep_frac * sigma_eff / sigma;
    
    %------------ photons ---------------- %
    % total number of photons available
    n_pho = num_src*num_pho;

    % number of photons per bayesian iteration
    n_pho_group = round(n_pho/50);
    n_pho_group = 2e4;
    
    % number of photons for direct detection pre-estimate
    n_imag_mu = 0;
    
    
    % ------------ Loop Through Trials ------------------ %
    
    % for each configuration run a certain number of reconstruction trials
    %parpool(num_workers)
    %parfor t = 1:trials
    for t = 1:trials
        
        % set seed
        seed = (array_id-1) * trials + t;
        
        % generate a random scene
        centroid_aligned = 1;
        src_brites = ones(num_src,1)/num_src;
        src_coords = genMinDistConstellation(src_brites, min_sep_frac, centroid_aligned);
        
        % scene [source brightnesses, source coordinates (in fractional rayleigh units)]
        scene = [src_brites, src_coords];
        
                
        % Run Adaptive Bayesian Multi-Parameter Estimation Algorithm
        [est_scene,est_trace] = measurement_multipho(scene, ...
                             'n_max',num_src, ...
                             'n_pho_group', n_pho_group, ...
                             'n_pho_SLD', n_pho, ...
                             'n_imag_mu', n_imag_mu,...
                             'n_pri',num_src,...
                             'n_sam',2e5,...
                             'bri_known',1,...
                             'proj_method','Personick',...
                             'per_eps', 0,...% 0<=per_eps<=1; HG 0 mode background
                             'n_HG_modes', max_order, ...
                             'aperture', aperture, ...
                             'U',U, ...
                             'dark_lambda', dark_pho_lambda, ...
                             'phase_sigma', phase_err_sigma, ...
                             'seed',seed...
                            ); 
        
                        
        % parameter estimates
        est_brites = est_scene(:,1);
        est_coords = est_scene(:,2:3);
              
        % compute localization error
        err = LocalizationError(src_coords, est_coords,1); 
        
        % data stucture for configuration
        cfg_data(t).n_pho = n_pho;
        cfg_data(t).n_pho_group = n_pho_group;
        cfg_data(t).scene = scene;
        cfg_data(t).est_scene = est_scene;
        cfg_data(t).est_trace = est_trace;
        cfg_data(t).err = err;

    end
    
    % save the configuration data to the global data structure
    DS.data(cfg_id{:}) = {cfg_data};
    
    % save the global data structure
    fname = [num2str(array_id),'_cfg','.mat'];
    save(fullfile(DS.save_dir,fname),'cfg_id','DS')
    

end