function xy = splitting(m, mu, sigma, min_xy, max_xy, varargin)

% Parser

   defaultmethod = 'ran_cen';
   defaultr_rat = 1.5;
   
   p = inputParser;

   addRequired(p,'m');
   addRequired(p,'mu');
   addRequired(p,'sigma');
   addRequired(p,'min_xy');
   addRequired(p,'max_xy');
   
   addOptional(p,'method',defaultmethod);
   addOptional(p,'r_rat',defaultr_rat);
   
   parse(p, m, mu, sigma, min_xy, max_xy, varargin{:});
   
   m = p.Results.m;
   mu = p.Results.mu;
   sigma = p.Results.sigma;
   min_xy = p.Results.min_xy;
   max_xy = p.Results.max_xy;
   
   method = p.Results.method;
   r_rat = p.Results.r_rat;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch method
    
    case 'ran_cen' % here sigma is a vector

while 0 == 0
    
        
     x(2:size(m,1),1) = normrnd(mu(1),sqrt(sigma(1,1))/8, size(m,1)-1,1);
     y(2:size(m,1),1) = normrnd(mu(2),sqrt(sigma(2,2))/8, size(m,1)-1,1);
     
     x(1) = ( sum(m)*mu(1) - sum( m.*x, 1 ) )/m(1);
     y(1) = ( sum(m)*mu(2) - sum( m.*y, 1 ) )/m(1);
     
     if and( and( x(1)>=min_xy(1), x(1)<=max_xy(1)) , ...
             and( y(1)>=min_xy(2), y(1)<=max_xy(2)) )
         
        break;
         
     else
         
         clear x y
         
         continue;
         
     end
    
    
end

    xy = [x,y];

case 'ran_cir' % here sigma is a matrix
        
    
     m = m/sum(m);
     
     [V,D] = eig(sigma);
     
     for i = 1:size(m,1)
        
         xy(1,i) = cos(2*pi*(i-1)/size(m,1));
         
         xy(2,i) = sin(2*pi*(i-1)/size(m,1));
         
     end
     
     xy = (D*V'*xy)'/2*r_rat + repmat( mu', [size(m,1), 1]); %/4 reduce variance
        
        
end


end
