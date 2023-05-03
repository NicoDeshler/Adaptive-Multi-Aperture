function measurement(varargin)

%%%%%%%%%%%%%%%%% Paser %%%%%%%%%%%%%%%%%%%%%%%%%%

    defaultseed = 0;
    defaultn_exp = 1;
    defaultn_imag_mu = 500;
    defaultn_pho_SLD = 9500;
    defaultscene = [];
    defaultuse_SLD = 1;
    defaultn_s = 8;
   
    p = inputParser;
    
    addOptional(p,'seed',defaultseed);
    addOptional(p,'n_exp',defaultn_exp);
    addOptional(p,'n_imag_mu',defaultn_imag_mu);
    addOptional(p,'n_pho_SLD',defaultn_pho_SLD);
    addOptional(p,'scene',defaultscene);
    addOptional(p,'use_SLD',defaultuse_SLD);
	addOptional(p,'n_s',defaultn_s);
   
    parse(p, varargin{:});
   
    seed = p.Results.seed;
    n_exp = p.Results.n_exp;
    n_imag_mu = p.Results.n_imag_mu;
    n_pho_SLD = p.Results.n_pho_SLD;
    scene = p.Results.scene;
    use_SLD = p.Results.use_SLD;
	n_s = p.Results.n_s;

%%%%%%%%%%%%%%%%%%%%%%% Scene generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global sig n_max

    sig = [1;1]; % doubled the noramlizated scale; don't change at the moment
	
	n_max = 6;
    
    if isempty(scene) 
    
        scene = scene_gen('seed', seed, ...
					      'n_s_mu', n_s, ...
						  's_random', 0);
        
    end
    
    if use_SLD == 1
    
        temp = imag_pho_gen(scene, 'seed', seed, 'n_exp', n_exp, 'n_imag_mu', n_imag_mu);
        
        pho_imag = temp{1};
    
        pho_imag_fil = temp{2};
    
        n_imag = temp{3};

        %%%%%%%%%%%%%%%%%%%%%%%% Imaging measurement
    
        EM_out = EM( n_imag, n_max, pho_imag, pho_imag_fil,'out',...
                    'EM = {}; EM{1} = [permute(m, [1,3,2]), tp]; EM{2} = r_ind''; ');
        tp = EM_out{1};
        pho_gp = EM_out{2};
    
    else
        
        temp = imag_pho_gen(scene, 'seed', seed, 'n_exp', n_exp, 'n_imag_mu', n_imag_mu+n_pho_SLD);
        
        pho_imag = temp{1};
    
        pho_imag_fil = temp{2};
    
        n_imag = temp{3};

        %%%%%%%%%%%%%%%%%%%%%%%% Imaging measurement
    
        imag_est = EM( n_imag, n_max, pho_imag, pho_imag_fil,'out',...
                      'EM = [permute(m, [1,3,2]), tp];');
                  
        for i = 1:n_exp
            
            imag_est(:,:,i) = sortrows(imag_est(:,:,i), 'descend');
            
        end
        
        clearvars -except imag_est seed scene tp
        
        filename = ['Imag_Est_',int2str(seed),'.mat'];

        save(filename);
        
        return;

    end
    



for i = 1:n_exp
    
    pho_gp_temp = zeros(1, size(pho_gp,2));
    
    tp_temp = tp(:,:,i);
    
    tp_temp = [(1:n_max)', tp_temp];
    
    tp_temp = sortrows(tp_temp, 2, 'descend');
    
    for j = 1:nnz(tp_temp(:,2))
       
        pho_gp_temp = pho_gp_temp + j*(pho_gp(i,:) == tp_temp(j,1));
        
    end
    
    pho_gp(i,:) = pho_gp_temp;
    
    tp(:,:,i) = tp_temp(:,2:4);
    
    clear tp_temp
    
end

%%%%%%% Proprotional representative system: grouping and splitting %%%%%%%%

% Even prior

rng(seed)

cand = {};

tic;
for ind = 1:n_exp
    
    cand_temp = tp(:,:,ind);
    
    cand_temp = cand_temp(1:nnz(tp(:,1,ind)),1:3);
    
    for j = size(cand_temp,1):n_max
        
        fprintf(['Here is grouping',int2str(j),'\n'])
        
        
        m_prior = ones(j,1)/j;
            
        temp = ini_splitting( tp(:,:,ind), pho_gp(ind, 1:n_imag(ind)),...
                              pho_imag(:, 1:n_imag(ind) ,ind), m_prior );
                              
        temp = zero_padding(temp);
            
        if j == size(cand_temp,1)
            
            cand{ind} = temp;
                
        else
                
            cand{ind} = cat(3, cand{ind}, temp);
                
        end
        
    end
end
toc;
%{
filename = ['SLD_Est_',int2str(seed),'.mat'];

save(filename);

return;
%}

%%%%%%%%%%%%%%%%%%%%%% SLD projector  %%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Eigenvector Generation

likelihood = zeros(n_max,n_exp);

mod_pick = zeros(n_max,6,n_exp);

mod_vote = zeros(n_max,n_exp);;

cand_check = {};

for ii = 1:n_pho_SLD
    
    clc;
    
    ii

prob ={};

for k = 1:n_exp
    
    if or( ii == 1, ~any(likelihood(1:nnz(likelihood(:,k)),k)))
        
        mod_pick(:,:,k) = model_sel(cand{k});
        
    else
        
        likelihood_temp = likelihood(1:nnz(likelihood(:,k)),k);
        
        vote = find( likelihood_temp == max(likelihood_temp));
        
        mod_vote(vote,k) = mod_vote(vote,k) + 1;
       
        mod_pick(:,:,k) = model_sel(cand{k},'method','Likelihood','likelihood',likelihood_temp);
        
    end
    
    tp_temp = mod_pick(:,:,k); 
    
    tp_temp = tp_temp(1:nnz(tp_temp(:,1)),:);
    
    n_now = size(tp_temp,1);
    
    eig_v = projectors(tp_temp);
    
    dtp_temp(:,1) = eig_v(1:n_now,1); 
    dtp_temp(:,2) = eig_v(n_now+1:2*n_now,1); 
    
    dtp(:,:,k) = zero_padding(dtp_temp);
     
    dm_temp = eig_v(2*n_now+1:3*n_now,1); 
     
    dm(:,k) = zero_padding(dm_temp);
    % Probability Generation
    prob_temp = prob_gen_3N(tp_temp(:,1), dm_temp, scene(:,1), ...
                            tp_temp(:,3:4), dtp_temp, scene(:,2:3));
    prob{k} = [max(1-sum(prob_temp), 0),prob_temp]; % 1-sum(prob) for no photon
    prob_sum(ii) = sum(prob_temp);
    
    % MIST
    pho_SLD(k) = pho_gen(prob{k}, 1);
    
    clear dtp_temp dm_temp
    
end


% Prior update
for k = 1:n_exp
    
    if pho_SLD(k) ~= 1
    
        temp = pri_up_3N(pho_SLD(k)-1, cand{k}, mod_pick(:,:,k), dtp(:,:,k), dm(:,k));
        
        cand{k} = temp{1};
        likelihood(:,k) = zero_padding(temp{2});
        
        if mod(ii,100) == 0
            
            cand_check{k,end+1} = cand{k};
            
        end

    else
        
        if mod(ii,100) == 0
            
            cand_check{k,end+1} = cand{k};
            
        end
        
        continue;
        
        
    end
    
    
    

end

clear L H eigQFI eig_v dtp dm mod_pick

end


% final model selection 

SLD_Est = {};

for k = 1:n_exp
    
    likelihood_temp = likelihood(1:nnz(likelihood(:,k)),k);
    
    SLD_Est{k} = model_final(cand{k},likelihood_temp);
    
    SLD_Est_vote{k} = model_final(cand{k},likelihood_temp,...
                                  'method', 'vote', ...
                                  'vote', mod_vote(:,k));

end

% clearvars -except cand SLD_Est seed scene tp

filename = ['SLD_Est_',int2str(seed),'.mat'];

save(filename);

end
