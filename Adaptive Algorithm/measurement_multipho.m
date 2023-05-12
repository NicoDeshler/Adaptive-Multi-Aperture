function [est, cand_check] = measurement_multipho(scene, varargin)

%%%%%%%%%%%%%%%%% Parser %%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultseed = 0;
    defaultn_imag_mu = 500;
    defaultn_pho_SLD = 9500;
    defaultuse_SLD = 1;
    defaultn_max = 10;
    defaultn_HG_modes = 5;
    defaultn_pho_group = 500;
    defaultn_pri = 0;
    defaultn_sam = 250000;
    defaultbri_known = 0;
    defaultproj_method = 'Personick';
    defaultper_eps = 0;
    defaultsavepath = 'out\est.mat'; % save filepath
    defaultaperture = [0,0,1];
    defaultU = 1;
    defaultdark_lambda = 0;
    defaultphase_sigma = 0;
    
    p = inputParser;
    
    addRequired(p,'scene');
    
    addOptional(p,'seed',defaultseed);
    addOptional(p,'n_imag_mu',defaultn_imag_mu);
    addOptional(p,'n_pho_SLD',defaultn_pho_SLD);
    addOptional(p,'use_SLD',defaultuse_SLD);
    addOptional(p,'n_max',defaultn_max);
    addOptional(p,'n_HG_modes',defaultn_HG_modes);
    addOptional(p,'n_pho_group',defaultn_pho_group);
    addOptional(p,'n_pri',defaultn_pri);
    addOptional(p,'n_sam',defaultn_sam);
    addOptional(p,'bri_known',defaultbri_known);
    addOptional(p,'proj_method',defaultproj_method);
    addOptional(p,'per_eps',defaultper_eps);
    addOptional(p,'savepath',defaultsavepath);
    addOptional(p,'aperture',defaultaperture);
    addOptional(p,'U',defaultU);
    addOptional(p,'dark_lambda',defaultdark_lambda);
    addOptional(p,'phase_sigma',defaultphase_sigma);
    
   
    parse(p, scene, varargin{:});
   
    scene = p.Results.scene;                % ground truth scene [Nx3] matrix
    seed = p.Results.seed;                  % random number seed
    n_imag_mu = p.Results.n_imag_mu;        % Poisson mean photon number for direct imaging 
    n_pho_SLD = p.Results.n_pho_SLD;        % Poisson mean photon number for adaptive imaging (total photon number)
    use_SLD = p.Results.use_SLD;            % SLD flag [0: Just direct imaging, 1: Adaptive Modal Imaging]  
    n_max = p.Results.n_max;                % max number of ources
    n_HG_modes = p.Results.n_HG_modes;      % number of HG modes
    n_pho_group = p.Results.n_pho_group;    % number of photons per adaptation
    n_pri = p.Results.n_pri;                % boolean flag indicating whether or not the number of sources is known apriori
    n_sam = p.Results.n_sam;                % number of posterior samples
    bri_known = p.Results.bri_known;        % 
    proj_method = p.Results.proj_method;    %
    per_eps = p.Results.per_eps;            % 
    savepath = p.Results.savepath;          % save filepath
    aperture = p.Results.aperture;          % multiaperture configuration
    U = p.Results.U;                        % mixing matrix for multiple apertures
    dark_lambda = p.Results.dark_lambda;    % Poisson mean of dark current photon noise
    phase_sigma = p.Results.phase_sigma;    % stdev of gaussian phasing error for multiple apertures 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n_pho_used = 0;
    
    if n_imag_mu > 0

    % Generate arrival photons on image plane
    
    pho_imag = imag_pho_gen(scene,...
                            'seed', seed,...
                            'n_imag_mu', n_imag_mu, ...
                            'per_eps', per_eps);
        

   %%%%%%%%%%%%%%%%%%%%%%%% Imaging measurement
    % Estimate the scene using direct imaging
    EM_out = EM( n_max, ...
                 pho_imag, ...
                 'out', 'EM = cell(1,2); EM{1} = [m, tp]; EM{2} = r_ind;');
    tp = EM_out{1};       % initial scene estimate          
    pho_gp = EM_out{2};   % number of sources associated with each measurement outcome
    
    if use_SLD == 0 % only imaging measurement
        
        %%%%%%%%%%%%%%%%%%%%%%%% Imaging measurement
    
        imag_est = sortrows(tp, 'descend');
        
        clearvars -except imag_est seed scene tp
        
        filename = ['Imag_Est_',int2str(seed),'.mat'];

        save(filename);
        
        return;

    end
    
    else
        
        % No imaging; initialization at (0,0)
        
        tp = zero_padding([1,0,0], n_max);
        pho_gp = 1;
        pho_imag = [0; 0];
        
    end

%%%%%%%%%%%%%% N_min = 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tp = zeros(n_max, 3);
tp(1,1) = 1;
tp(1,2:3) = mean(pho_imag,2);
pho_gp = ones(1,size(pho_imag,2));

%%%%%%% Proprotional representative system: grouping and splitting %%%%%%%%

% Even prior

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(seed)

tic;

cand_temp = tp(1:nnz(tp(:,1)),:);

n_model = n_max-size(cand_temp,1)+1;

cand = zeros(n_max, 6, n_model);%%%%%%%%%%%%%%%%%%%%% chenge here %%%%%%%%%%%%%%%%%%

for j = size(cand_temp,1):n_max %%%%%%%%%%%%%%%%%%%%% chenge here %%%%%%%%%%%%%%%%%%
        
        fprintf(['Here is grouping ',int2str(j),'\n'])
        
        m_prior = ones(j,1)/j;
            
        temp = ini_splitting( tp, ...
                              pho_gp, ...
                              pho_imag, ...
                              m_prior, ...
                              'gp_method', 'simple_even');
                          
        posCov{j-size(cand_temp,1)+1} = diag([temp(:,5);temp(:,6)])/8;
            
        cand(:,:,j-size(cand_temp,1)+1) = zero_padding(temp, n_max); %%%%%%%%%%%%%%%%%%%%% chenge here %%%%%%%%%%%%%%%%%%
        
end


%%%%%%%%%%%%%%%%%%% prior_info %%%%%%%%%%%%%%%%%%%%%%%%

cand = prior_info(cand, scene, ...
                 'n_pri', n_pri, ...
                 'pos_pri', 0, ...
                 'bri_pri', 0, ...
                 'pos_var_pri', 0, ...
                 'bri_var_pri', 10);
 

%_---------------------
%%%%% NICO EDIT  %%%%%%
%_---------------------
cand(:,3:4) = normrnd(0,1e-3,n_pri,2); % new initialization of positions
cand(:,5:6) = 1e-2; % new initialization of positions

    

if n_pri
    
    posCov_temp = posCov;
    posCov{1} = posCov_temp{end}; % need to put n_s = n_max

end

%%%%%%%%%%%%%%%%%%% prior_info %%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%% SLD projection  %%%%%%%%%%%%%%%%%%%%

n_SLD_true = poissrnd(n_pho_SLD);

likelihood = zeros(n_model,1);
likelihood_check = [];

mod_vote = zeros(n_model,1);

cand_check{1} = cand; 

%%%%%%%%%%%%%%%%%%%%%%%% Above checked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ii = 1; %period
% This is the bayesian update loop. Keep running until all photons used
while n_pho_used < n_SLD_true
    
    ii
    
    pho_now = poissrnd(n_pho_group);
    
    pho_now = min(pho_now, n_SLD_true - n_pho_used);
    
    
%%%%%%%%%%%%%%%%%%%%%% Eigenvector Generation %%%%%%%%%%%%%
    
    if ii == 1 || ~any(likelihood)
        
        mod_pick = model_sel(cand);
        
    else
        
        vote = find( likelihood == max(likelihood));
        
        mod_vote(vote) = mod_vote(vote) + 1;
       
        mod_pick = model_sel(cand, ...
                            'method','Likelihood',...
                            'likelihood',likelihood);
        
    end
    
    tp= mod_pick(1:nnz(mod_pick(:,1)),:);
    
    L = SLD_delta_HG(n_HG_modes, tp, ...
                     'method', proj_method, ..., ...
                     'bri_flag', 0,...
                     'aperture',aperture,...
                     'U',U...
                     );
    
    % Probability Generation
    prob = prob_gen_2nd_mom(scene, L, ...
                            'n_modes', n_HG_modes,...
                            'aperture', aperture,...
                            'U',U .* exp(1i*normrnd(0,phase_sigma,1,size(U,2)))); % add piston phasing error to mixing matrix;
                        
    
    % MIST: modify the number of the photon
    pho_SLD = pho_gen(prob,  pho_now);
    
    % sample dark current photons
    pho_dark = pho_gen_dark(size(L,1),dark_lambda);
    
    % collect all measured photons
    pho_SLD = [pho_SLD, pho_dark];

    
% Prior update

    if ~isempty(pho_SLD)
        
        temp = Bay_2nd_mom(pho_SLD, cand, L, n_HG_modes, n_sam , n_max , ii,...
                          'posCov', [],...
                          'bri_known', bri_known, ... 
                          'cv_rate', [1,1,50],...
                          'aperture',aperture,...
                          'U',U);
                
                
        
        % Centred at original coordinates
        cand = temp{1};

        cand_check{end+1} = cand;
        
        likelihood = temp{2};
            
        posCov = temp{3};
        
        likelihood_check(:,ii) = likelihood;
        
    end
    
    clear dtp dm

    ii = ii + 1;
    n_pho_used = n_pho_used + pho_now;
    
    toc;
    
    %save(savepath,'scene','cand_check','likelihood_check')
    
end

toc;

% final model selection 
est = model_final(cand,likelihood);
%save(savepath,'scene','est','cand_check','likelihood_check')

end
