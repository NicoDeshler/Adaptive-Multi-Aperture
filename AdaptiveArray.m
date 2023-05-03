function AdaptiveArray(cfg_id)



DS = DSAdaptive();

trials = 50;        % trials per configuration
max_order = 5;
min_sep_frac = 2.^linspace(-6,-3,4);
num_src = [3,4,5,6];
flux = [1e3,1e4,1e5];
apertures = {};
aperture_names = {};
rl = [];

[a,n,m,p] = cfg_id = 




% generate a random scene

for t = 1:DS.trials

    % generate a random scene
    centroid_aligned = 0;
    scene_brites = ones(num_src,1)/num_src;
    scene_coords_frac = genMinDistConstellation(src_brites, min_sep, centroid_aligned);
    
    
    
    
    
    % Run Adaptive Bayesian Multi-Parameter Estimation Algorithm
    scene_est  = measurement_multipho(scene, ...
                         'n_max',num_src, ...
                         'n_pho_group', 10000, ...
                         'n_pho_SLD', 500000, ... % number of photons used throughout bayesian update
                         'n_imag_mu', 1000,...
                         'n_pri',num_src,...
                         'bri_known',1,...
                         'proj_method','Personick',...
                         'per_eps', 0,...% 0<=per_eps<=1; HG 0 mode background
                         'n_HG_modes', max_order, ...
                         'aperture', aperture,...
                         'U',U...
                        ); 
     scene_est 
                
end
                
end