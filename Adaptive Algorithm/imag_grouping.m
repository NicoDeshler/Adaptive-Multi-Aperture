function imag_grouping(I_gp, count, history, leftover, m_prior, varargin)

   defaultgp_method = 'full';

   p = inputParser;
   
   addRequired(p,'I_gp');
   addRequired(p,'count');
   addRequired(p,'history');
   addRequired(p,'leftover');
   addRequired(p,'m_prior');
   
   addOptional(p,'gp_method',defaultgp_method);
   
   parse(p, I_gp, count, history, leftover, m_prior, varargin{:});
   
   I_gp = p.Results.I_gp;
   count = p.Results.count;
   history= p.Results.history;
   leftover= p.Results.leftover;
   m_prior = p.Results.m_prior;
   gp_method = p.Results.gp_method;
   

global group del_gp;


switch gp_method
    
    case 'simple_even'
        
    m_temp = m_prior(1);
    
    I_gp_temp = I_gp;
    
    for i = 1:size(I_gp,1)
       
        while I_gp_temp(i) >= m_temp 
            
            group = [group; i];
            
            I_gp_temp(i) = I_gp_temp(i) - m_temp;
            
        end
        
    end
    
    if size(group,1) < size(m_prior,1)
        
        I_gp_temp = [I_gp_temp, (1:size(I_gp_temp,1))'];
        
        I_gp_temp = sortrows(I_gp_temp,'descend');
        
        group = [group; I_gp_temp(1:(size(m_prior,1) - size(group,1)),2)];
        
    end
    
    case 'full' 

   m_temp = repmat(m_prior(count),[size(I_gp,1),1]);
   
   gp = (I_gp - m_temp >= -del_gp).*(I_gp ~= 0);
   
   gp = find(gp);

if size(m_prior,1) ~= count
    
    if ~isempty(gp)
    
        for i = 1:size(gp,1)
       
            [count i];
       
            history(count) = gp(i);
            
            I_gp_temp = I_gp;
       
            I_gp_temp(gp(i)) = (I_gp_temp(gp(i)) - m_prior(count));
       
            I_gp_temp(gp(i)) = I_gp_temp(gp(i))*(I_gp_temp(gp(i)) >= 0);
            
            count_temp = count + 1;
       
            imag_grouping(I_gp_temp,count_temp,history,leftover,m_prior);
       
        end
        
    else
        
            leftover = [leftover; count];
        
            count_temp = count + 1;
       
            imag_grouping(I_gp,count_temp,history,leftover,m_prior);     
        
        
    end
    
else
    
    switch sum(del_gp,1)
        
        case 0
            
            if ~isempty(gp)
    
                for i = 1:size(gp,1)
                    
                    history_temp = history;
                    
                    I_gp_temp = I_gp;
                    
                    history_temp(count) = gp(i);
       
                    I_gp_temp(gp(i)) = (I_gp_temp(gp(i)) - m_prior(count));
                    
                    I_gp_temp = sortrows([(1:size(I_gp_temp,1))', I_gp_temp] ,2, 'descend');
                    
                    for j = 1:size(leftover,1)
                    
                        history_temp(leftover(j)) = I_gp_temp(j,1);
                    
                    end
                    
                    group = [group, history_temp];
                
                end
                
        
            else
                
                leftover = [leftover; count];
                    
                I_gp = sortrows([(1:size(I_gp,1))', I_gp] ,2, 'descend');
                
                for j = 1:size(leftover,1)
                    
                    history(leftover(j)) = I_gp(j,1);
                    
                end
        
                group = [group, history];
        
            end
        
        otherwise
            
            if ~isempty(gp)
    
                for i = 1:size(gp,1)
            
                    I_gp(gp(i)) = I_gp(gp(i)) - m_prior(count);
            
                    history(count) = gp(i);
                    
                    group = [group, history];
                
                end
        
            else
        
                group = [group, history];
        
            end

    end
    
end

end

end