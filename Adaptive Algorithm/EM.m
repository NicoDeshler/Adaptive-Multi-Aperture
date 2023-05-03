function EM = EM(n_max, pho, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   defaultN_em = 1000;
   defaultOut = strcat('EM = [m, tp];');
   defaultsig = [1,1];
   defaultrange_x = max(defaultsig)*4;
   defaultfix_m = [];

   p = inputParser;

   addRequired(p,'n_max');
   addRequired(p,'pho');

   addOptional(p,'n_em',defaultN_em);
   addOptional(p,'out',defaultOut);
   addOptional(p,'sig',defaultsig);
   addOptional(p,'range_x',defaultrange_x);
   addOptional(p,'fix_m',defaultfix_m);
   
   
   parse(p, n_max, pho, varargin{:});
   
   n_max = p.Results.n_max; % just the average photon
   pho = p.Results.pho;
   
   n_em = p.Results.n_em;
   out = p.Results.out;
   sig = p.Results.sig;
   range_x = p.Results.range_x;
   fix_m = p.Results.fix_m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% initializaton %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(-range_x,range_x,1001);
dx = x(2) - x(1);
[X,Y] = meshgrid(x,x);

n_imag = size(pho,2);

pho_mean = sum(pho,2)./size(pho,2);

pho_var = sum( (pho - repmat(pho_mean, [1,n_imag]) ).^2, 2)./size(pho,2);

rat_var = 1;

% even initialization
m = ones(n_max,1)/n_max;


    for i = 1:n_max
        % row: points; column: x,y%    
        tp(i,1) = pho_var(1)*cos(2*pi*(i-1)/n_max) + pho_mean(1);
        
        tp(i,2) = pho_var(2)*sin(2*pi*(i-1)/n_max) + pho_mean(2);

    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;

tp_old = zeros(size(tp));
m_old = zeros(size(m));

%while count<n_em
while ( ~isempty(find(tp-tp_old)) || ~isempty(find(m-m_old))) && count<n_em
    
%count
    
tp_old = tp;
m_old = m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% E-step

%%%%%%%%%%%%%%%% Change below

% dim_1: x,y; dim_2: photons; dim_3: point sources
    
    tp_temp = permute( repmat( tp, [1,1,n_imag] ), [2,3,1]);
    
    pho_temp = repmat( pho, [1,1,size(tp,1)] );

    r = permute( repmat(m, [1,1,n_imag] ), [2,3,1])...
       .*( exp( -(pho_temp(1,:,:)-tp_temp(1,:,:)).^2 ...
                -(pho_temp(2,:,:)-tp_temp(2,:,:)).^2  )...
        * sqrt(sig(1)*sig(2)*2/pi) ).^2*dx^2;      


[~,r_ind] = max(r,[],3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% M-step

for i = 1:n_max % i:point ; k:iteration 
    
    m(i) = sum( r_ind == i, 2 );
    
    tp(i,:,:) = [sum( (r_ind == i).*pho(1,:), 2)/m(i),...
                 sum( (r_ind == i).*pho(2,:), 2)/m(i)];  
                
end

if isempty(fix_m)
    
    m = m/n_imag;
    
else

    m = fix_m;
    
end

count = count + 1;

tp(isnan(tp)) = 0;

end

if count == n_em
    disp('max_n_em');
end

eval(out);


end