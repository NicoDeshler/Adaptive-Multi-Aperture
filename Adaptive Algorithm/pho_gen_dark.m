function pho = pho_gen_dark(n_modes,lambda)
    % generate dark current photons counts with rate lambda    
    pho = [];
    
    for k = 1:n_modes+1
        pho = [pho, repmat(k,1,poissrnd(lambda))];
    end
        
end