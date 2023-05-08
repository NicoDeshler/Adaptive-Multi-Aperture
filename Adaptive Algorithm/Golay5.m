function golay5 = Golay5(R_max)
    th = linspace(0,2*pi,5);
    r = [0,R_max*ones(1,4)];
    
    [x,y] = pol2cart(th,r); 
    golay5 = [x',y'];
end
