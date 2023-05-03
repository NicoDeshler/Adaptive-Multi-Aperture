function pho_gen = pho_gen(prob, n_pho)    

addpath('../');

n = size(prob,2); 

n_exp = size(prob,1); 

prob = prob/sum(prob);

prob = [zeros(n_exp,1), prob];
    
CDF = zeros(n_exp,n+1);
    
temp = rand([n_exp, n_pho]);
    
samp = zeros([n_exp, n_pho]);

for i = 1:n
    
   CDF(:,i+1) = sum(prob(:,1:i+1),2);
   samp_temp1 = repmat(CDF(:,i),[1, n_pho])< temp;
   samp_temp2 = temp <= repmat(CDF(:,i+1),[1, n_pho]);
   samp_temp = samp_temp1.*samp_temp2;
   
   samp = samp + i*samp_temp;
end

pho_gen = samp;

end
