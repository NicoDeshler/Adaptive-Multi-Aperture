function golay7= Golay7(R_max)
R1 = R_max * ones(1,3);
R2 = R_max/2 * ones(1,3);

th1 = pi/2 + linspace(0,2/3*2*pi,3);
th2 = th1 + (2*pi/3)/3;

[x,y] = pol2cart([th1,th2],[R1,R2]);
golay7 = [[0;x'],[0;y']];

    
end