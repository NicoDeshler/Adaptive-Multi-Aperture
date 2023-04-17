% dim is the dimension of the problem
% m is a dim X 1 vetor, with the relative intensities of the sources
% tp is the position of the sources, with diemnsion n X dim
% sig is the width of PSF; not used for dim = 1

function SLD = SLD_3N(m,tp,dim,varargin)


   defaultsig = [1;1];
   defaultOut = strcat('SLD = L;');
   defaultAp = [0,0];
   defaultdtp = [0,0];
   defaultdm = [0,0];
   defaultI_norm = 1;

   p = inputParser;
   
   addRequired(p,'m');
   addRequired(p,'tp');
   addRequired(p,'dim');
   
   addOptional(p,'dtp',defaultdtp);
   addOptional(p,'dm',defaultdtp);
   addOptional(p,'sig',defaultsig);
   addOptional(p,'out',defaultOut);
   addOptional(p,'Ap',defaultAp);
   addOptional(p,'I_norm',defaultI_norm);
   
   parse(p,m,tp,dim,varargin{:});
   
   m = p.Results.m;
   tp = p.Results.tp;
   dtp= p.Results.dtp;
   dm= p.Results.dm;
   dim = p.Results.dim;
   sig = p.Results.sig;
   out = p.Results.out;
   ap = p.Results.Ap;
   I_norm = p.Results.I_norm;
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if dim == 20
    
% Dimension 2

% Number of sources
n = size(m,1);
x = tp(:,1);
y = tp(:,2);

[X1,X2] = meshgrid(x,x);
[Y1,Y2] = meshgrid(y,y);


% Method 1
% Gram matrix
G = prod(sig)*exp(-(X1-X2).^2/2-(Y1-Y2).^2/2);

% rho in original basis
rho = G*diag(m)*G;

[V_rho, D_rho] = eig(rho,G);

[V_rho, D_rho] = sortem(V_rho,D_rho); 

% Subspace orthogonalization
% Group the subspaces
Subspace = {};
del = 1e-12; % can be set to other small number
d_rho = diag(D_rho);

count = 1;
temp = [1];

for i = 1:n-1

   if abs(d_rho(i+1)-d_rho(i))<del
       
      temp = [temp,i+1];
      
      if i == n-1
           Subspace{count} = temp;
      end
      
   else
       
       Subspace{count} = temp;
       count =  count + 1;
       clear temp;
       temp = [i+1];
       
       if i == n-1
           Subspace{count} = temp;
       end
   end
    
end

V_num = 0;

for i = 1:numel(Subspace)
    
    if numel(Subspace{i}) == 1
        
        continue;
        
    else
        
        V_rho(:,Subspace{i}) = GramSchmidt(V_rho(:,Subspace{i}));
        
    end
   
end

% Normalization in original basis
for i = 1:n
   %pk(i) = sum(sum(diag(V_rho(:,i))*G*diag(V_rho(:,i)),1),2);
   V_rho(:,i) = V_rho(:,i)/sqrt( sum(sum(diag(V_rho(:,i))*G*diag(V_rho(:,i)),1),2) );
    
end
%{

% Method 2
% Gram matrix
G = prod(sig)*exp(-(X1-X2).^2/2-(Y1-Y2).^2/2);

[V_G, D_G] = eig(G);

% rho in original basis
rho = G*diag(m)*G;

% transformation of rho
t_rho = sqrt(D_G)*V_G'*rho*V_G*sqrt(D_G);

% eigvalues and eigvectors of t_rho
[V_rho, D_rho] = eig(t_rho);

% tramsformation of V_rho
V_rho = D_G*V_G'*V_rho;
%}

% Normalization in original basis
for i = 1:n
   %pk(i) = sum(sum(diag(V_rho(:,i))*G*diag(V_rho(:,i)),1),2);
   V_rho(:,i) = V_rho(:,i)/sqrt( sum(sum(diag(V_rho(:,i))*G*diag(V_rho(:,i)),1),2) );
    
end

% from here
% dG and Aper

dG = zeros(n,n,2*n);
Aper = zeros(n,n,2*n);

for k = 1:n % derivatives; In genera should be N X M; N #sources M #parameters
    
    for i = 1:n
        for j = 1:n
       
            if k == j
                
               %dG(i,j,k) = 0;
               dG(i,j,k) = (x(i)-x(j))*G(i,j); % check i - j or j - i
               dG(i,j,n+k) = (y(i)-y(j))*G(i,j); % check i - j or j - i 
                
            end
        end
    end
    Aper(:,:,k) = inv(V_rho)*inv(G)*dG(:,:,k);
    Aper(:,:,n+k) = inv(V_rho)*inv(G)*dG(:,:,n+k);
end

% d2G and Apara

d2G_4D = zeros(n,n,2*n,2*n);
d2G = {};
Aper_cb = {};

for k = 1:2*n % #parameters
    for l = 1:2*n % #parameters
        
        for i = 1:n
            for j = 1:n
                
                if and(mod(k,n) == mod(i,n), mod(l,n) == mod(j,n))
                %if and(k == i, l == j)
                    if and(k<=n, l<=n)
                    
                        d2G_4D(i,j,k,l) = (1-(x(i)-x(j))^2)*G(i,j);
                        
                    elseif and(k>n, l>n)
                        
                        d2G_4D(i,j,k,l) = (1-(y(i)-y(j))^2)*G(i,j);
                        
                    else
                        
                        d2G_4D(i,j,k,l) = -(x(i)-x(j))*(y(i)-y(j))*G(i,j);
                
                    end
                    
                end
                
            end
        end
        d2G{k,l} = d2G_4D(:,:,k,l);
        Aper_cb{k,l} = (Aper(:,:,k))'*Aper(:,:,l);
    end
end



Apara = cell2mat(d2G) - cell2mat(Aper_cb);
[U,S,V] = svd(Apara);
Apara = U*sqrt(S)*V';
% Further reduce the space; No need in general
%{
columnsWithAllZeros = all(abs(Apara')/max(max(abs(Apara))) < 1e-14);
Apara = Apara(~columnsWithAllZeros, :);
%}


temp = zeros(2*n^2,n,2*n);
for i = 1:2*n
   
    temp(:,:,i) = Apara(:,(i-1)*n+1:i*n);
    
end

clear Apara

Apara = temp;

% SLD 

L = zeros(2*n^2+n,2*n^2+n,3*n);  % Not adaptive

[D1, D2] = meshgrid(diag(D_rho),diag(D_rho));

for k = 1:2*n % #parameters
   
    L(1:n,1:n,k) = 2*( Aper(:,:,k)*diag(m)*G*V_rho ...
                    + (Aper(:,:,k)*diag(m)*G*V_rho)' ) ...
                   ./(D1+D2);
     
    L(n+1:2*n^2+n,1:n,k) ...
        = (2*Apara(:,:,k)*diag(m)*G*V_rho)./repmat((diag(D_rho))',[2*n^2,1]);
    
    L(1:n,n+1:2*n^2+n,k) = L(n+1:2*n^2+n,1:n,k)';
     
end

GC = G*V_rho;

for k = 1:n % #parameters
    
    v_GC = GC(k,:);
   
    L(1:n,1:n,2*n+k) = 2*( v_GC'*v_GC ...
                         - GC'*diag(m)*GC ) ...
                       ./(D1+D2);
     
end

L(1:n,1:n,3*n+1) = D_rho;
% to here

eval(out);

%SLD = L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif dim == 21
    
% Dimension 2

% Number of sources
n = size(m,1);

% Number of Aperture
n_ap = size(ap,1);

% Normalization of coordinates

x = tp(:,1);
y = tp(:,2);

%{
x = tp(:,1);
y = tp(:,2);
%}
dx = squeeze(dtp(:,1,:));
dy = squeeze(dtp(:,2,:));

[X2,X1] = meshgrid(x,x);
[Y2,Y1] = meshgrid(y,y);

% Method 1
% Gram matrix
G_temp = prod(sig)*exp(-(X1-X2).^2/2-(Y1-Y2).^2/2);
G = zeros(size(X1));

for i = 1:size(ap,1)
    
    G = G + G_temp.*exp(1i*ap(i,1)*(X1-X2)).*exp(1i*ap(i,2)*(Y1-Y2))/n_ap;
    
end

% rho in original basis
rho = G*diag(m)*G/I_norm;

[V_rho, D_rho] = eig(rho,G);

[V_rho, D_rho] = sortem(V_rho,D_rho); 

%V_rho = V_rho./ repmat(sqrt(sum(V_rho.*conj(V_rho),1)),[n,1]);

V_rho_test = V_rho; 
% Subspace orthogonalization

% Group the subspaces
Subspace = {};
del = 1e-12; % can be set to other small number
d_rho = diag(D_rho);

count = 1;
temp = [1];

for i = 1:n-1

   if abs(d_rho(i+1)-d_rho(i))<del
       
      temp = [temp,i+1];
      
      if i == n-1
           Subspace{count} = temp;
      end
      
   else
       
       Subspace{count} = temp;
       count =  count + 1;
       clear temp;
       temp = [i+1];
       
       if i == n-1
           Subspace{count} = temp;
       end
   end
    
end

V_num = 0;

for i = 1:numel(Subspace)
    
    if numel(Subspace{i}) == 1
        
        continue;
        
    else
        
        V_rho(:,Subspace{i}) = GramSchmidt(V_rho(:,Subspace{i}));
        
    end
   
end

% Normalization in original basis
for i = 1:n
   
   pk(i) = sum(sum(diag(V_rho(:,i))*G*diag(V_rho(:,i)),1),2);
   V_rho(:,i) = V_rho(:,i)/sqrt( real(sum(sum(diag(V_rho(:,i))*G*diag(V_rho(:,i)),1),2)) );
    
end

% dG and Aper

dG = zeros(n,n,size(dtp,3));
Aper = zeros(n,n,size(dtp,3));

for k = 1:size(dtp,3) % derivatives; In general should be N X M; N #sources M #parameters
    
    for i = 1:n
        for j = 1:n
       
            for l = 1:size(ap,1)
            
            dG(i,j,k) = dG(i,j,k) ...
                      + G_temp(i,j)...
                      .*exp(1i*ap(l,1)*(x(i) - x(j)))...
                      .*exp(1i*ap(l,2)*(y(i) - y(j)))...
                      .*( dx(j,k)*(x(i) - x(j) - 1i*ap(l,1)) ...
                        + dy(j,k)*(y(i) - y(j) - 1i*ap(l,2)) )/n_ap;      
            
            end
        end
    end
    
    Aper(:,:,k) = inv(V_rho)*inv(G)*dG(:,:,k);

end

% d2G and Apara

d2G_4D = zeros(n,n,size(dtp,3),size(dtp,3));

d2G = {};
Aper_cb = {};

for k = 1:size(dtp,3) % #parameters
    for l = 1:size(dtp,3) % #parameters
        
        for i = 1:n
            for j = 1:n
                
                for p = 1:size(ap,1)

            d2G_4D(i,j,k,l) = d2G_4D(i,j,k,l) ...
                            + G_temp(i,j)...
                            .*exp(1i*ap(p,1)*(x(i) - x(j)))...
                            .*exp(1i*ap(p,2)*(y(i) - y(j)))...
                            .*( ... 
                              [dx(i,k), dy(i,k)]*( eye(2) ...
                            - [x(i)-x(j);y(i)-y(j)]*[x(i)-x(j),y(i)-y(j)] )...
                            * [dx(j,l); dy(j,l)] ...
                            ...
                            + ( ap(p,1)*dx(i,k) + ap(p,2)*dy(i,k) ) ...
                             *( ap(p,1)*dx(j,l) + ap(p,2)*dy(j,l) ) ...
                            ...
                            +1i*(...
                                ( ap(p,1)*dx(i,k) + ap(p,2)*dy(i,k) )...
                                * ( dx(j,l)*(x(i)-x(j)) + dy(j,l)*(y(i)-y(j)) )...
                               +( ap(p,1)*dx(j,l) + ap(p,2)*dy(j,l) )...
                                * ( dx(i,k)*(x(i)-x(j)) + dy(i,k)*(y(i)-y(j)) )...
                                ) ...
                            )/n_ap; 
                end
                
            end
        end
        d2G{k,l} = d2G_4D(:,:,k,l);
        Aper_cb{k,l} = (Aper(:,:,k))'*Aper(:,:,l);
    end
end


del = 1e-24;
%{
[V,S] = eig(cell2mat(d2G) - cell2mat(Aper_cb));
diag(S)

Apara = V*sqrt( S.*(S>=0) + (S<0)*del)*V';
eig(Apara)

Apara = sqrtm(cell2mat(d2G) - cell2mat(Aper_cb));

eig(cell2mat(d2G) - cell2mat(Aper_cb))
%}
Apara = vpa(cell2mat(d2G) - cell2mat(Aper_cb));
[U,S,V] = svd(Apara);
Apara = vpa(U*sqrt(S)*V');


% Further reduce the space; No need in general
%{
columnsWithAllZeros = all(abs(Apara')/max(max(abs(Apara))) < 1e-14);
Apara = Apara(~columnsWithAllZeros, :);
%}

s_sub = size(dtp,3)*n;

temp = zeros(s_sub,n);

for i = 1:size(dtp,3)
   
    temp(:,:,i) = Apara(:,(i-1)*n+1:i*n);
    
end

clear Apara

Apara = temp;

% SLD 

L = zeros(s_sub+n,s_sub+n,size(dtp,3));  % Not adaptive

[D1, D2] = meshgrid(diag(D_rho),diag(D_rho));

GC = G*V_rho;

for k = 1:size(dtp,3) % #parameters
    
    L(1:n,1:n,k) = 2*( ...
                       GC'*( diag(dm(:,k)) - sum(dm(:,k),1)*diag(m) )*GC ...
                    +  Aper(:,:,k)*diag(m)*G*V_rho ...
                    + (Aper(:,:,k)*diag(m)*G*V_rho)' ) ...
                   ./(D1+D2);
     
    L(n+1:s_sub+n,1:n,k) ...
        = (2*Apara(:,:,k)*diag(m)*G*V_rho)./repmat((diag(D_rho))',[s_sub,1]);
    
    L(1:n,n+1:s_sub+n,k) = L(n+1:s_sub+n,1:n,k)';
     
end

L(1:n,1:n,size(dtp,3)+1) = D_rho;
% to here

eval(out);

%SLD = L;
end