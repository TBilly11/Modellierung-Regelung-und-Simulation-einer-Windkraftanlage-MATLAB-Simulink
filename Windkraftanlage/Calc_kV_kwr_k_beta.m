%num = (1:1:20)';
% Define the constants and parameters
c1  = 0.005;
c2	= 1.53;
c3	= 0.5;
c4	= 0.18;
c5	= 121;
c6	= 27.9;
c7	= 198;
c8	= 2.36;
c9	= 5.74;
c10	= 11.35;
c11	= 16.1;
c12	= 201;
rho = 1.225; % air density in kg/m^3
R = 63; % rotor radius in m
J = para_mdl.Jr + para_mdl.ngear^2 * para_mdl.Jg;
%k_Opt = 0.5*para_mdl.rho*pi*(para_mdl.R^5)*(para_mdl.Cp_max/((para_mdl.ngear^3)*(para_mdl.lambda_opt^3))); % Define k_Opt here

% New values for beta and v
beta = [0.0000, 0.0306, 0.0560, 0.0772, 0.0953, 0.1113, 0.1254, 0.1382, 0.1499, 0.1607, 0.1783, 0.2134, 0.2421, 0.2669, 0.2890, 0.3091, 0.3279, 0.3456, 0.3625, 0.3789]';
v = [11.2600, 11.5378, 11.8156, 12.0933, 12.3711, 12.6489, 12.9267, 13.2044, 13.4822, 13.7600, 14.2600, 15.4533, 16.6467, 17.8400, 19.0333, 20.2267, 21.4200, 22.6133, 23.8067, 25.0000]';

omega_r = 1.2671; % rated rotor speed in rad/s

% Define the auxiliary functions
lambda_i = @(lambda,beta) 1./(lambda+0.08*beta)-0.035./(c11+c12*beta.^3);
f1 = @(lambda) c4./lambda;
f2 = @(lambda,beta) (c5*lambda_i(lambda,beta))-c6*beta-c7*beta.^c8-c9;
f3 = @(lambda,beta) exp(-c10*lambda_i(lambda,beta));

% Define the analytical CQ function and its derivatives
dlambda_i_dlambda = @(lambda,beta) -1./(lambda+0.08*beta).^2;

cQ_tilde = @(lambda,beta) c1*(1+c2*(beta+c3).^2)+f1(lambda).*f2(lambda,beta).*f3(lambda,beta);
cQ = @(lambda,beta) cQ_tilde(lambda,beta).*(1+sign(cQ_tilde(lambda,beta)))./2;
dcQ_dlambda = @(lambda,beta) -c4./lambda.^2.*f2(lambda,beta).*f3(lambda,beta)+f1(lambda).*(c5.*dlambda_i_dlambda(lambda,beta)).*f3(lambda,beta)+f1(lambda).*f2(lambda,beta).*(-c10.*exp(-c10*lambda_i(lambda,beta)).*dlambda_i_dlambda(lambda,beta));
dcQ_dbeta = @(lambda,beta) f1(lambda).*(c5.*dlambda_i_dbeta(lambda,beta)-c6-c7*c8*beta.^(c8-1)).*f3(lambda,beta)+f1(lambda).*f2(lambda,beta).*(-c10.*exp(-c10*lambda_i(lambda,beta)).*dlambda_i_dbeta(lambda,beta))+0.5*c1*c2*(beta+c3).^(-0.5);

% Define the derivatives of lambda_i function
dlambda_i_dlambda = @(lambda,beta) -1./(lambda+0.08*beta).^2;
dlambda_i_dbeta = @(lambda,beta) -0.08./(lambda+0.08*beta).^2+(3*0.035*c12*beta.^2)./(c11+c12*beta.^3).^2;

% Define the lambda function as a function of wind speed and rotor speed
lambda = @(V,omega_r) omega_r*R./V;

% Calculate the three formulas for each wind speed
k_omega = zeros(length(v),1); % initialize the vector for k_omega
k_beta = zeros(length(v),1); % initialize the vector for k_beta
k_V = zeros(length(v),1); % initialize the vector for k_V
for i = 1:length(v)
    k_omega(i) = 0.5*rho*v(i)^2*pi*R^3*(R/v(i)*dcQ_dlambda(lambda(v(i),omega_r),beta(i))); % calculate k_omega for the i-th wind speed
    k_beta(i) = 0.5*rho*v(i)^2*pi*R^3*dcQ_dbeta(lambda(v(i),omega_r),beta(i)); % calculate k_beta for the i-th wind speed
    k_V(i) = 0.5*rho*pi*R^3*(2*v(i)*cQ(lambda(v(i),omega_r),beta(i))-omega_r*R*dcQ_dlambda(lambda(v(i),omega_r),beta(i))); % calculate k_V for the i-th wind speed
end

% Create a table with the results
%T = table(v, k_omega, k_beta, k_V, 'VariableNames', {'WindSpeed', 'k_omega', 'k_beta', 'k_V'});

% Write the table to an Excel file
%writetable(T, 'output.xlsx');
% Define the numerator and denominator coefficients
% num = [(-8.6/J) 0]; % s
% den = [0 1 -(-2.63/J)]; % s^2 + 3s + 2

% Create the transfer function model
% G = tf (num,den)
tau_ref = 3.5;
a = -k_omega/J;
b = k_beta/J;
k_P = 1./(b.*tau_ref);
k_I = a./(b.*tau_ref);
%-----------------------------------------

index_kP = find(k_P == out.KP1.Data);
index_ki = find(k_I == out.KI1.Data);
% index_kbeta = find(k_beta == out.k_beta1.Data);
% index_komega = find(k_omega == out.k_omega1.Data);

fprintf('Index in k_P: %d\n', index_kP);
fprintf('Index in k_i: %d\n', index_ki);
% fprintf('Index in k_beta: %d\n', index_kbeta);
% fprintf('Index in k_omega: %d\n', index_komega);
% % Create a table with the results
% T = table(v, k_omega, k_beta, k_V,k_P, k_I,  'VariableNames', {'WindSpeed', 'k_omega', 'k_beta', 'k_V','k_P','k_I'});
% 
% %Write the table to an Excel file
% writetable(T, 'output1.xlsx');
