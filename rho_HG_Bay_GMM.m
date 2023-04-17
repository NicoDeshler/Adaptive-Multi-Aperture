function rho_HG = rho_HG_Bay_GMM(n_modes, scene, V, Mom, varargin)


   % Parser
   defaultbri_flag = 0;
   
   p = inputParser;

   addRequired(p,'n_modes');
   addRequired(p,'scene');
   addRequired(p,'V');
   addRequired(p,'Mom');
   
   addOptional(p,'bri_flag',defaultbri_flag );
   
   parse(p, n_modes, scene, V, Mom, varargin{:});
   
   n_modes = p.Results.n_modes;
   scene = p.Results.scene;
   V = p.Results.V;
   Mom = p.Results.Mom;
   
   bri_flag = p.Results.bri_flag;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    scene = sortrows(scene, 'descend');

    scene = scene(1:nnz(scene(:,1)), :);
    
    %N_modes = 3*n_modes*(3*n_modes+1)/2;
    N_modes = n_modes;
    
    if Mom == 1
        
        if bri_flag
            h = [V(1:size(scene,1)), V(size(scene,1)+1:2*size(scene,1))];
            b = V( 2*size(scene,1)+1:3*size(scene,1) );
            
        else
            h = [V(1:size(scene,1)), V(size(scene,1)+1:2*size(scene,1))];
        end
            
    end
    
    alp0 = sum(scene(:,2),1);
    alp = scene(:,2)/alp0;
    Alp = diag(alp)/(alp0+1) + alp*alp'*( 1 - 1/(alp0+1) );

    load 'aperture.mat'
    n_ap = size(aperture,1);
    H(:,:,1) = aperture(:,1)-aperture(:,1)';
    H(:,:,2) = aperture(:,2)-aperture(:,2)';
    
    num_sources =  size(scene,1);

    % instantiate the density operator or personick operators
    rho_HG = zeros(N_modes,N_modes);
    
    % posterior parameters
    bar_mu  = (2 * psf_sig) * scene(:,3:4);              % mean of gaussian posteriors for each source
    bar_sig = (2 * psf_sig) * scene(:,5:6);             % standard deviation of gaussian posteriors for each source
    %bar_mu  = psf_sig * scene(:,3:4);
    %bar_sig = psf_sig * scene(:,5:6);
    s_b = scene(:,1);                   % source brightness estimates
    
    % GAMMA0
    if Mom == 0
        for j1 = 1:N_modes         % linear HG union aperture index 
            for j2 = 1:N_modes     % linear HG union aperture index 

                p1 = pj(j1);
                q1 = qj(j1);
                nu1 = uj(j1);

                p2 = pj(j2);
                q2 = qj(j2);
                nu2 = uj(j2);

                n1 = p1 + p2;
                n2 = q1 + q2;

                for k = 1:size(scene,1)
                    for mu1 = 1:n_ap
                        for mu2 = 1:n_ap
                            b = H(mu1,mu2,:); % baseline components
                            rho_HG(j1,j2) = rho_HG(j1,j2) + 1/n_ap * s_b(k) * conj(U(nu1,mu1))*U(nu2,mu2) * ...
                                                            Chi(n1,p1,p2,b(1),psf_sig,bar_mu(k,1),bar_sig(k,1)) *...
                                                            Chi(n2,q1,q2,b(2),psf_sig,bar_mu(k,2),bar_sig(k,2));
                        end
                    end
                end
            end
        end
        % trace for Gamma 0 should equal 1
    
    % GAMMA 1
    elseif Mom == 1
        
        for j1 = 1:N_modes         % linear HG union aperture index 
            for j2 = 1:N_modes    % linear HG union aperture index 

                p1 = pj(j1);
                q1 = qj(j1);
                nu1 = uj(j1);

                p2 = pj(j2);
                q2 = qj(j2);
                nu2 = uj(j2);

                n1 = p1 + p2;
                n2 = q1 + q2;

                for k = 1:size(scene,1)

                    % moment sum
                    h_mu = sum(h(1:num_sources~=k,:).*bar_mu(1:num_sources~=k,:),'all');

                    % summations
                    for mu1 = 1:n_ap
                        for mu2 = 1:n_ap
                            b = H(mu1,mu2,:); % baseline components
                            rho_HG(j1,j2) = rho_HG(j1,j2) + 1/n_ap * s_b(k) * conj(U(nu1,mu1))*U(nu2,mu2) * (...
                                                            h(k,1)*Chi(n1+1,p1,p2,b(1),psf_sig,bar_mu(k,1),bar_sig(k,1)) * ...
                                                            h(k,2)*Chi(n2+1,q1,q2,b(2),psf_sig,bar_mu(k,2),bar_sig(k,2)) + h_mu * ...
                                                            Chi(n1,p1,p2,b(1),psf_sig,bar_mu(k,1),bar_sig(k,1)) * ...
                                                            Chi(n2,q1,q2,b(2),psf_sig,bar_mu(k,2),bar_sig(k,2)) ...
                                                            );
                        end
                    end
                end
            end
        end
        % trace for Gamma 1 should equal sum(h.*bar_mu,'all')
        
        
    end

end


% helper function
function chi_element = Chi(n,p1,p2,b,psf_sig,bar_mu,bar_sig)
% computes the closed form intermediate integrals of a multi-gaussian psf and a
% gaussian posterior.
% Inputs:
% n - the gaussian moment 
% p1 - a 1D HG index
% p2 - a 1D HG index
% b  - a componenent of a baseline
% psf_sig - the gaussian PSf standard deviation for an individual sub-aperture
% bar_mu - the mean of the gaussian posterior
% bar_sig - the standard deviation of the gaussian posterior


%{
b = b * 2*psf_sig;
sig_2 = 2*bar_sig.^2 + 1;
a = sqrt(bar_sig.^2 ./ sig_2);
t = (bar_mu - 1i*b.*bar_sig.^2)/sig_2;
A = 1./sqrt(factorial(p1) .* factorial(p2) .* 2*pi*bar_sig.^2) .* exp( - 0.5 * ( ( bar_mu - 1i* b ).^2 ./  sig_2  + b.^2 ) );



%psf_sig = 1/sqrt(2);
psf_sig = 1;
%bar_sig = bar_sig./psf_sig;

%sig_2 = bar_sig.^2 + 2*psf_sig.^2;
sig_2 = 2*bar_sig.^2 + psf_sig.^2;

A = 1./sqrt(factorial(p1) .* factorial(p2) .* 4*pi*bar_sig.^2) .* exp( - 0.5 * ( ( bar_mu - 1i* b ).^2 ./  sig_2  + b.^2 ) );
a = sqrt( 4 * bar_sig.^2  ./ sig_2 );
t = (bar_mu - 4*1i*b*bar_sig.^2) ./ sig_2;

%}


sig_2 = bar_sig.^2 + 2*psf_sig.^2;

A = 1./sqrt(factorial(p1).*factorial(p2)) .* 1/sqrt(2*pi*bar_sig.^2) .* exp( - ( ( bar_mu - 2*1i* b.*psf_sig.^2 ).^2 ./ (2 * sig_2 ) + b.^2 .* psf_sig.^2 ) );
a = sqrt( 2*bar_sig.^2 .* psf_sig.^2 ./ sig_2 );
t = 2*psf_sig.^2 .* (bar_mu - 1i*b*bar_sig.^2) ./ sig_2;

moments = compute_gaussian_moments(t, a, n);
x_mom_n = moments(n+1);

chi_element = A.*a.*sqrt(2*pi).*x_mom_n;

%chi_element = A .* sqrt(2*pi) .* a .* (1i*a).^n .* exp(- (t ./ (2*a) ).^2) .* pu(n, -1i * t ./ a);
end
