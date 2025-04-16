%----------------------------------------------------------------------------------
% plot_Cx_BEM.m
%----------------------------------------------------------------------------------
% Description: Plot the CT, CQ and CP look up table over lambda as tip speed ratio 
%              and with the parameter beta as pitch angle 
%-------------------------------------------------------------------------------------

%load para_mdl_WT.mat             
%load Cx_maps.mat



figure

% plot CQ
subplot(221)
plot(lambda_array,CQ(:,:))
title('Torque coefficients:  $C_Q(\lambda,\beta)$','Interpreter', 'latex')
xl = xlabel('$\lambda$', 'Interpreter', 'latex');
yl = ylabel('$C_Q$', 'Interpreter', 'latex');
set(xl, 'Fontsize', 14);
set(yl, 'Fontsize', 16);
grid

% plot CT   
subplot(222)
plot(lambda_array,CT(:,:))
title('Thrust coefficients: $C_T(\lambda,\beta)$','Interpreter', 'latex')
xl = xlabel('$\lambda$', 'Interpreter', 'latex');
yl = ylabel('$C_T$', 'Interpreter', 'latex');
set(xl, 'Fontsize', 14);
set(yl, 'Fontsize', 16);
grid     

% plot CP
subplot(223)
plot(lambda_array,CP(:,:))
title('Power coefficients: $C_P(\lambda,\beta)$','Interpreter', 'latex')
xl = xlabel('$\lambda$', 'Interpreter', 'latex');
yl = ylabel('$C_P$', 'Interpreter', 'latex');
set(xl, 'Fontsize', 14);
set(yl, 'Fontsize', 16);
grid


% plot rotor blade chord 
% figure(4)
% plot(BlSec.r,BlSec.t,'-o')
% title('blade chord','Interpreter', 'latex')
% xl = xlabel('Blade span [m]', 'Interpreter', 'latex');
% yl = ylabel('Chord [m]', 'Interpreter', 'latex');
% set(xl, 'Fontsize', 14);
% set(yl, 'Fontsize', 16);
% grid
% 
% 
% % plot rotor blade chord 
% figure(5)
% plot(BlSec.r,BlSec.alpha_twist,'-o')
% title('blade twist angle','Interpreter', 'latex')
% xl = xlabel('Blade span [m]', 'Interpreter', 'latex');
% yl = ylabel('alpha twist', 'Interpreter', 'latex');
% set(xl, 'Fontsize', 14);
% set(yl, 'Fontsize', 16);
% grid

