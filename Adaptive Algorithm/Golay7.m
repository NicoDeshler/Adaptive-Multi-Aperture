function golay7= Golay7(R)
% reference:
% Manx Arrays: Perfect Non-Redundant Interferometric Geometries
% McKay et. al. 2022
% doi:10.1029/2022RS007500


e1 = [1;0];
e2 = [cos(pi/3);sin(pi/3)];

M = [0,0;
     2,0;
     1,1;
     0,-2;
     1,-2;
     -2,1;
     -2,2];
 
golay7 = M*[e1,e2]';

% shift array to be centroid aligned
golay7 = golay7 - mean(golay7);

% rescale to have R_max radius
golay7 = golay7 * R ./ max(vecnorm(golay7,2,2));

end