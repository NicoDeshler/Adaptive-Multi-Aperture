function moments = compute_gaussian_moments(mu, sigma, max_order)
% Uses a recurrence relation to compute raw moments of the gaussian
% distribution.
% 
% mu: mean of the Gaussian distribution
% sigma: standard deviation of the Gaussian distribution
% max_order: the maximum order of the moments to compute

% Thank you GPT4!

% Initialize the moments vector
moments = zeros(max_order+1,1);

% Compute the 0th order moment
moments(1) = 1;
% Compute the 1st order moment
moments(2) = mu;

% Compute the higher-order moments using the recurrence relation
for m = 2:max_order
    moments(m+1) = mu*moments(m) + (m-1)*sigma^2*moments(m-1);
end
end