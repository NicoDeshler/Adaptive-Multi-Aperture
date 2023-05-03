function cand = ini_splitting(tp, pho_gp, pho, m_prior, varargin)

   defaultgp_method = 'full';

   p = inputParser;
   
   addRequired(p,'tp');
   addRequired(p,'pho_gp');
   addRequired(p,'pho');
   addRequired(p,'m_prior');
   
   addOptional(p,'gp_method',defaultgp_method);
   
   parse(p, tp, pho_gp, pho, m_prior, varargin{:});
   
   tp = p.Results.tp;
   pho_gp = p.Results.pho_gp;
   pho = p.Results.pho;
   m_prior = p.Results.m_prior;
   gp_method = p.Results.gp_method;

global group del_gp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Grouping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    I_gp = tp(1:nnz(tp(:,1)), 1);
    
    del_gp = 0*I_gp;

    group = [];

    history = zeros(size(m_prior,1),1);

    count = 1;

    imag_grouping(I_gp, count, history, [], m_prior, 'gp_method', gp_method)

    % Deleteing degenergcy
    
    m_gp =  findgroups(m_prior);

    for i = 1:max(m_gp)
    
    ind = find(m_gp == i);
    
        if numel(ind) == 1
                continue;
        else
                group(ind,:) = sort(group(ind,:));
        end
    
    end 

    group = unique(group.','rows').';

%%%%%%%%%%%%%%%%%%%%%%% Set up candmeter matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cand is the candmeters: 
    % c1: intensity
    % c2: candmeters of Dirichlet distribution
    % c3: x mean
    % c4: y mean
    % c5: x variance
    % c6: y variance
    cand = zeros(size(m_prior,1), 6, size(group,2));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Splitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    for k = 1:size(group, 2) % possible situations
        
        m_norm = sum( m_prior.*(group(:,k)~=0), 1); % Group index
        
        cand_temp = [];
    
        for j = unique(nonzeros(group(:,k)))'
       
            m_temp = m_prior( group(:,k) == j )*I_gp(j) ...
                     /sum(m_prior( group(:,k) == j )) ...
                     /sum(I_gp(unique(nonzeros(group(:,k)))')); % Assign prior into each group
        
            %pho_temp = repmat(pho_gp == j, [2,1]).*pho; % coordinates of group photons
            
            %pho_temp( :, ~any(pho_temp,1) ) = [];
            
            pho_temp = pho(:, pho_gp == j);
            
            mu = mean(pho_temp,2);
            
            sigma2 = cov(pho_temp(1,:),pho_temp(2,:));
            
            if any(isnan(sigma2(:))) || rank(sigma2) ~=2
               
                sigma2 = eye(2)/4;
                
            end

            %%%%%%%%%%%% splitting the point %%%%%%%%%%%%%%
            
            on9(:,1) = m_temp;
                
            on9(:,2) = 1+m_temp;
            
            if numel(m_temp) == 1
                
                on9(3:4) = mu';
                
                on9(5:6) = (diag(sigma2))';
                
            else
                
                max(pho, [], 2);
                
                min(pho, [], 2);
                
                on9(:,3:4) = splitting(m_temp, mu, sigma2,...
                                       min(pho, [], 2), max(pho, [], 2), ...
                                       'method','ran_cir');
                
                on9(:,5:6) = repmat([sigma2(1,1),sigma2(2,2)],[size(m_temp,1),1]);
                
            end
            
            cand_temp = [cand_temp; on9];
            
            clear on9
        
        end
        
        cand(1:size(cand_temp,1),:,k) = cand_temp;
        
        
    
    end


end