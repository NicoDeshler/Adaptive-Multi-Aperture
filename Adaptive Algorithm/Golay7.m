function golay7= Golay7(R_max)
% reference:
% Manx Arrays: Perfect Non-Redundant Interferometric Geometries
% McKay et. al. 2022
% doi:10.1029/2022RS007500

% unit vectors in triangular coordinate grid
e1 = [1;0];
e2 = [cos(pi/3);sin(pi/3)];

M = [2,0;
     1,1;
     1,-1;
     0,-1;
     -1,1;
     -1,2];
 
golay7 = M*[e1,e2];

% shift array to be centroid aligned
golay7 = golay7 - mean(golay7);

% add center point
golay7 = [[0,0];golay7];

% rescale to have R_max radius
golay7 = golay7 * R_max ./ max(vecnorm(golay7,2));

%{
R1 = R_max * ones(1,3);
R2 = R_max/2 * ones(1,3);

th1 = pi/2 + linspace(0,2/3*2*pi,3);
th2 = th1 + (2*pi/3)/3;

[x,y] = pol2cart([th1,th2],[R1,R2]);
golay7 = [[0;x'],[0;y']];
%}

    
end