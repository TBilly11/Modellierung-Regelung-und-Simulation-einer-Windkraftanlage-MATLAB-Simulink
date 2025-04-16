%----------------------------------------------------------------------------------
% plot_chord_twist.m
%----------------------------------------------------------------------------------
% Description: Plot chord and twist of each blade elemment
%-------------------------------------------------------------------------------------

load BladeSec_NREL_5MW_Data.mat  % blade element data


BlSec.r           = BladeSec_Data(:,1);   % blade element distance to rotation axis
BlSec.alpha_twist = BladeSec_Data(:,2);   % aero twist [grad]
BlSec.dr          = BladeSec_Data(:,3);   % blade element length
BlSec.t           = BladeSec_Data(:,4);   % blade chord
BlSec.NFoil       = BladeSec_Data(:,5);   % airfoil ID (number) 


figure

% plot rotor blade chord 
subplot(121)
plot(BlSec.r,BlSec.t,'-o')
title('blade chord','Interpreter', 'latex')
xl = xlabel('Blade span [m]', 'Interpreter', 'latex');
yl = ylabel('Chord [m]', 'Interpreter', 'latex');
set(xl, 'Fontsize', 14);
set(yl, 'Fontsize', 16);
grid


% plot rotor twist 
subplot(122)
plot(BlSec.r,BlSec.alpha_twist,'-o')
title('blade twist angle','Interpreter', 'latex')
xl = xlabel('Blade span [m]', 'Interpreter', 'latex');
yl = ylabel('alpha twist', 'Interpreter', 'latex');
set(xl, 'Fontsize', 14);
set(yl, 'Fontsize', 16);
grid

