function B = phase_fn(xy_coords,aperture,U)
    % computes the phase function representing local mode intereference of
    % multiple apertures in a tensor product basis.
    Phi = exp(1i*aperture(:,1:2)*xy_coords);
    B = sqrt(1/size(U,1)) * U * Phi;
end