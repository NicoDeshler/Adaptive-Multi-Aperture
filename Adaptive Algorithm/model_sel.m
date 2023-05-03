function model = model_sel(tp, varargin)

    defaultMethod = 'random';
    defaultLikelihood = 0;
   
    p = inputParser;
   
    addRequired(p,'tp');
   
    addOptional(p,'method',defaultMethod);
    addOptional(p,'likelihood',defaultLikelihood);
   
    parse(p, tp, varargin{:});
   
    tp = p.Results.tp;
   
    method = p.Results.method;
    likelihood = p.Results.likelihood;
    
    
    switch method
        
        case 'random'
            
            model = tp(:,:,randi([1 size(tp,3)]));
            
        case 'Likelihood'
            
            [~,ind] = max(likelihood);
            
            model = tp(:,:,ind);
            
    end


end