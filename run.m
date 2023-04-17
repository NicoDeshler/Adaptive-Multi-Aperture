function run(varargin)

clc

%%%%%%%%%%%%%%%%% Paser %%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultn_case = 10;
    defaultn_imag_mu = 500;
    defaultn_pho_SLD = 9500;
    defaultn_exp = 1;
    defaultscene = [];
    defaultuse_SLD = 1;
   
    p = inputParser;
    
    addOptional(p,'n_case',defaultn_case);
    addOptional(p,'n_imag_mu', defaultn_imag_mu);
    addOptional(p,'n_pho_SLD', defaultn_pho_SLD);
    addOptional(p,'n_exp', defaultn_exp);
    addOptional(p,'scene',defaultscene);
    addOptional(p,'use_SLD',defaultuse_SLD);
   
    parse(p, varargin{:});
   
    n_case = p.Results.n_case; 
    n_imag_mu = p.Results.n_imag_mu;
    n_pho_SLD = p.Results.n_pho_SLD;
    scene = p.Results.scene;
    n_exp = p.Results.n_exp;
    use_SLD = p.Results.use_SLD;
    
    if ~isempty(scene)
        
        measurement('n_pho_SLD',n_pho_SLD,...
                    'n_exp',n_exp,...
                    'n_imag_mu',n_imag_mu,...
                    'scene',scene,...
                    'use_SLD',use_SLD);
        
    else
    
    parfor i = n_case
    
        measurement('seed',i,...
                    'n_pho_SLD',n_pho_SLD,...
                    'n_imag_mu',n_imag_mu,...
                    'n_exp',n_exp,...
                    'use_SLD',use_SLD);
    
    end
    
    end
    
    

end