load 'aperture.mat'
n_ap = size(aperture,1);
H(:,:,1) = aperture(:,1)-aperture(:,1)';
H(:,:,2) = aperture(:,2)-aperture(:,2)';
h;       % 2xN_sources projection vector;


% instantiate the density operator or its moments
rho_HG = zeros(N_modes,N_modes);

% GAMMA0
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
        
        s_b = scene(:,2);
        
        for k = 1:size(scene,1)
            
            bar_mu = scene(k,3:4);          % mean of gaussian posteriors for source k
            bar_sig = scene(k,5:6);         % standard deviation of gaussian posteriors for source k
            
            for mu1 = 1:n_ap
                for mu2 = 1:n_ap
                    b = H(mu1,mu2,:); % baseline components
                    rho_HG(j1,j2) = rho_HG(j1,j2) + 1/num_sources * s_b(k) * conj(U(nu1,mu1))*U(nu2,mu2) * ...
                                                    Chi(n1,p1,p2,b(1),psf_sig,bar_mu(k,1),bar_sig(k,1)) *...
                                                    Chi(n2,q1,q2,b(2),psf_sig,bar_mu(k,2),bar_sig(k,2));
                end
            end
        end
    end
end

% GAMMA 1
% GAMMA0
for j1 = 1:n_ap*N_modes*(N_modes+1)/2         % linear HG union aperture index 
    for j2 = 1:n_ap*N_modes*(N_modes+1)/2     % linear HG union aperture index 
        
        p1 = indices(1,j1);
        q1 = indices(2,j1);
        nu1 = indices(3,j1);
        
        p2 = indices(1,j2);
        q2 = indices(2,j2);
        nu2 = indices(3,j2);
        
        n1 = p1 + p2;
        n2 = q1 + q2;
        
        s_b = scene(:,2);
        
        for k = 1:size(scene,1)
            
            % moment sum
            h_mu = sum(h(i~=k,:).*bar_mu(i~=k,:),'all');
            
            % summations
            for mu1 = 1:n_ap
                for mu2 = 1:n_ap
                    b = H(mu1,mu2,:); % baseline components
                    rho_HG(j1,j2) = rho_HG(j1,j2) + 1/num_sources * s_b(k) * conj(U(nu1,mu1))*U(nu2,mu2) * (...
                                                    Chi(n1+1,p1,p2,b(1),psf_sig,bar_mu(k,1),bar_sig(k,1)) * ...
                                                    Chi(n2+1,q1,q2,b(2),psf_sig,bar_mu(k,2),bar_sig(k,2)) + h_mu * ...
                                                    Chi(n1,p1,p2,b(1),psf_sig,bar_mu(k,1),bar_sig(k,1)) * ...
                                                    Chi(n2,q1,q2,b(2),psf_sig,bar_mu(k,2),bar_sig(k,2)) ...
                                                    );
                end
            end
        end
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



sig_2 = bar_sig.^2 + 2*psf_sig.^2;

A = 1./sqrt(factorial(p1).*factorial(p2)) .* 1/sqrt(2*pi*bar_sig.^2) .* exp( - ( ( bar_mu - 2*1i* b.*sig_psf.^2 ).^2 ./ (2 * sig_2 ) + b.^2 .* psf_sig.^2 ) );
a = sqrt( 2*bar_sig.^2 .* psf_sig.^2 ./ sig_2 );
t = 2*psf_sig.^2 .* (bar_mu - 1i*b*bar_sig.^2) ./ sig_2;


chi_element = A .* sqrt(2*pi.*a) .* (1i*a).^n .* exp(- (t ./ (2*a) ).^2) .* pu(-1i * t ./ a);
end



               