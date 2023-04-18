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
    
    load 'aperture.mat'
    n_ap = size(aperture,1);
    H(:,:,1) = aperture(:,1)-aperture(:,1)';
    H(:,:,2) = aperture(:,2)-aperture(:,2)';
    
    [pj,qj,uj] = Indices_HG_GMMAperture(n_modes,n_ap);
    
    %N_modes = 3*n_modes*(3*n_modes+1)/2;
    N_modes = numel(pj);

    
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
    
    num_sources =  size(scene,1);

    % instantiate the density operator or personick operators
    rho_HG = zeros(N_modes,N_modes);
    
    % posterior parameters
    %bar_mu  = (2 * psf_sig) * scene(:,3:4);              % mean of gaussian posteriors for each source
    %bar_sig = (2 * psf_sig) * scene(:,5:6);             % standard deviation of gaussian posteriors for each source
    %bar_mu  = psf_sig * scene(:,3:4);
    %bar_sig = psf_sig * scene(:,5:6);
    bar_mu  = scene(:,3:4);
    bar_sig = scene(:,5:6);
    
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
                            rho_HG(j1,j2) = rho_HG(j1,j2) + 1/n_ap * s_b(k) * conj(U(nu1,mu1))*U(nu2,mu2)...
                                                            *Chi(n1,p1,p2,b(1),bar_mu(k,1),bar_sig(k,1))...
                                                            *Chi(n2,q1,q2,b(2),bar_mu(k,2),bar_sig(k,2));
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
                            rho_HG(j1,j2) = rho_HG(j1,j2) + 1/n_ap * s_b(k) * conj(U(nu1,mu1))*U(nu2,mu2)...
                                                            * (h(k,1)*Chi(n1+1,p1,p2,b(1),bar_mu(k,1),bar_sig(k,1))...
                                                            *  Chi(n2,q1,q2,b(2),bar_mu(k,2),bar_sig(k,2))...
                                                            +  h(k,2)*Chi(n2+1,q1,q2,b(2),bar_mu(k,2),bar_sig(k,2)) ...
                                                            *  Chi(n1,p1,p2,b(1),bar_mu(k,1),bar_sig(k,1))... 
                                                            + h_mu...
                                                            * Chi(n1,p1,p2,b(1),bar_mu(k,1),bar_sig(k,1))...
                                                            * Chi(n2,q1,q2,b(2),bar_mu(k,2),bar_sig(k,2)));
                        end
                    end
                end
            end
        end
        % trace for Gamma 1 should equal sum(h.*bar_mu,'all')
        
    end
end


% helper function
function chi_element = Chi(n,p1,p2,b,bar_mu,bar_sig)
%function chi_element = Chi(n,p1,p2,b,psf_sig,bar_mu,bar_sig)
% computes the closed form intermediate integrals of a multi-gaussian psf and a
% gaussian posterior.
% Inputs:
% n - the gaussian moment 
% b  - a componenent of a baseline
% psf_sig - the gaussian PSf standard deviation for an individual sub-aperture
% bar_mu - the mean of the gaussian posterior
% bar_sig - the standard deviation of the gaussian posterior

% primary variables from completing the square
xi = sqrt(2*bar_sig^2 + 1);
a = bar_sig/xi;                         % effective stdev
t = (bar_mu - 1i*b*bar_sig^2) / (xi^2); % effective mean

% secondary variables from completing the square
u = bar_mu/bar_sig;                     % second var
tau = -1i * b / (2 * bar_sig);          % second var mu
bar_a = sqrt(-xi^2 / (2 * bar_sig^2));  % second var std

% gaussian moment
moments = compute_gaussian_moments(t, a, n);
x_mom_n = moments(n+1);

% output element
chi_element = 1/xi * exp( 1/2 * ( ( (u-tau) / bar_a )^2 + (bar_a * a * b)^2 )) * x_mom_n; % this line outputs the same as Kwan's mom_Bay function
chi_element = chi_element / sqrt(factorial(p1)*factorial(p2));
end
