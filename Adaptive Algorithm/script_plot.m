clear all
clc

ind = [51,52,53,54];

ind_exp = 1;

for kk = ind
    
    
    str_load_Imag = ['load("Imag_Est_',int2str(kk),'.mat")'];
    str_load_SLD = ['load("SLD_Est_',int2str(kk),'.mat")'];
    
    eval(str_load_Imag);
    eval(str_load_SLD);

    % Display density and sampling on image plane

range_x = 4;
x = linspace(-range_x,range_x,1001);
dx = x(2) - x(1);
[X,Y] = meshgrid(x,x);

m = scene(:,1);

tp0 = scene(:,2:3);

I = zeros(size(X));

for i = 1:size(tp0,1)
    
   I = I + m(i)*( exp( -(X-tp0(i,1)).^2 -(Y-tp0(i,2)).^2  )...
               *sqrt(2/pi) ).^2*dx^2;  

end

figure             
surf(X,Y,I,'EdgeColor','none')
xlabel('x');
ylabel('y');
title('Focal Plane Intensity')
colormap jet;
colorbar;
view(0,90)
set(gca,'FontSize',14);
axis tight;

rel_size = 1000;

% Initial Guess
tp_temp = tp;

tp_temp( ~any(tp_temp,2), : ) = [];

% Imag

imag_temp = imag_est(:,:,ind_exp);

imag_temp = imag_temp(1:nnz(imag_temp(:,1)),:);


% SLD
SLD_temp = SLD_Est{ind_exp};

SLD_temp = SLD_temp(1:nnz(SLD_temp(:,1)),:);

% plot
figure('units','normalized','outerposition',[0 0 1 1])
   
scatter( scene(:,2),scene(:,3),scene(:,1)*rel_size,'k','filled')
hold on;

scatter( tp_temp(:,2),tp_temp(:,3),tp_temp(:,1)*rel_size,'y','filled')
hold on;

scatter( imag_temp(:,2),imag_temp(:,3),imag_temp(:,1)*rel_size,'b','LineWidth',2)
hold on;

scatter( SLD_temp(:,2),SLD_temp(:,3),SLD_temp(:,1)*rel_size,'rs','LineWidth',2)
hold on;
   

xlabel('x');
ylabel('y');
axis([-max(x), max(x), -max(x), max(x)]);
set(gca,'FontSize',30);
legend('Ground Truth', 'Initial Guess', 'Imaging','SLD Proj', 'Location', 'northwest')

    
    
end