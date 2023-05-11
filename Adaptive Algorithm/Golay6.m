function golay6 = Golay6(R)
% reference:
% Manx Arrays: Perfect Non-Redundant Interferometric Geometries
% McKay et. al. 2022
% doi:10.1029/2022RS007500

% unit vectors in triangular coordinate grid
e1 = [1;0];
e2 = [cos(pi/3);sin(pi/3)];

M = [2,0;
     1,1;
     0,-2;
     1,-2;
     -2,1;
     -2,2];
 
golay6 = M*[e1,e2]';

% shift array to be centroid aligned
golay6 = golay6 - mean(golay6);

% rescale to have R_max radius
golay6 = golay6 * R ./ max(vecnorm(golay6,2,2));

end