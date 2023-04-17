 function prob_gen = prob_gen_3N(m, dm, m0, tp, dtp, tp0, varargin)

format short
 
% Parser

   defaultDim = 21;
   
   p = inputParser;

   addRequired(p,'m');
   addRequired(p,'dm');
   addRequired(p,'m0');
   addRequired(p,'tp');
   addRequired(p,'dtp');
   addRequired(p,'tp0');
   
   addOptional(p,'dim',defaultDim);
   
   parse(p, m, dm, m0, tp, dtp, tp0, varargin{:});
   
   m = p.Results.m;
   dm = p.Results.dm;
   m0 = p.Results.m0;
   tp = p.Results.tp;
   dtp = p.Results.dtp;
   tp0 = p.Results.tp0;
   
   dim = p.Results.dim;
   
%addpath('../');

n = size(m,1);
n0 = size(m0,1);

L = SLD_3N(m,tp,dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = L;');
C = SLD_3N(m,tp,dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = V_rho;');
Aper = SLD_3N(m,tp,dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = Aper;'); % Chnage for multi parameters
Apara = SLD_3N(m,tp,dim, 'dm', dm, 'dtp', dtp, 'out', 'SLD = Apara;');

Omega = -Aper/Apara;

[V_L, D_L] = eig(L(:,:,1)); % Chnage for multi parameter
[V_L, D_L] = sortem(V_L, D_L);

l_rho = V_L(1:n,:);
l_mu = V_L(n+1:end,:);

%{
% for normalization check
G = SLD(m,tp,dim,'dtp',dtp, 'out', 'SLD = G;');
dG = SLD(m,tp,dim,'dtp',dtp, 'out', 'SLD = dG;');
dG2 = SLD(m,tp,dim,'dtp',dtp, 'out', 'SLD = cell2mat(d2G);');

%Normalization checking

N_check = {};

% dot product in rho space; should be n X n identity 
N_check{1,1} = C'*G*C; 

% dot product in rho and mu space; should be 0
N_check{1,2} = C'*dG/Apara + C'*G*C*Omega; 

% dot product in rho and mu space; should be 0
N_check{2,1} = N_check{1,2}'; 

% dot product in mu space; should be n X n identity 

N_check{2,2} =  Apara'\dG2/Apara + Omega'*C'*dG/Apara ...
              + (Omega'*C'*dG/Apara)' + Omega'*C'*G*C*Omega;


N_check = cell2mat(N_check)


on91 = Apara'\dG2/Apara
on92 = Omega'*C'*dG/Apara
on93 = (Omega'*C'*dG/Apara)'
on94 = Omega'*C'*G*C*Omega

on99 = N_check{2,2}



%eig(inv(Apara)*Apara);

%eig(N_check);

%} 

% For distance (\gamma)
[X, X0] = meshgrid(tp(:,1),tp0(:,1));
[Y, Y0] = meshgrid(tp(:,2),tp0(:,2));

% G = prod(sig)*exp(-(X1-X2).^2/2-(Y1-Y2).^2/2);
G0 = exp(-(X0-X).^2/2-(Y0-Y).^2/2);

for i = 1:n0
        for j = 1:n
            
            dG0(i,j) = G0(i,j)...
                    .*( dtp(j,1)*(tp0(i,1) - tp(j,1)) ...
                      + dtp(j,2)*(tp0(i,2) - tp(j,2)) );  
        end
end

prob_gen = zeros(1,size(L,2));

%{
Roh = {};

Rho{1,1} = (G0*C)'*diag(m0)*G0*C;
Rho{1,2} = (G0*C)'*diag(m0)*( dG0/Apara + G0*C*Omega );
Rho{2,1} = Rho{1,2}';
Rho{2,2} = ( dG0/Apara + G0*C*Omega )'*diag(m0)*( dG0/Apara + G0*C*Omega );

rho = cell2mat(Rho);

eig(rho);
%}

for i = 1:size(L,2)
        
        prob_gen(i) = m0'*abs( G0*C*l_rho(:,i) ...
                     + ( dG0/Apara + G0*C*Omega)*l_mu(:,i) ).^2;
          %{       
        prob_gen(i) = (V_L(:,i))'*rho*V_L(:,i);
        %}
end


end