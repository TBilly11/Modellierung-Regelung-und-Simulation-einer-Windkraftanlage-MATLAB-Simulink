%clear all

load para_mdl_WT_with_aero.mat

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

J = para_mdl.Jr + para_mdl.ngear^2 * para_mdl.Jg;
% New values for beta and v


omega_r = para_mdl.omega_r_R;

beta   = [0.0000,  0.0306,  0.0560,  0.0772,  0.0953,  0.1113,  0.1254,  0.1382,  0.1499,  0.1607,  0.1783,  0.2134,  0.2421,  0.2669,  0.2890,  0.3091,  0.3279,  0.3456,  0.3625,  0.3789]';
v      = [11.2600, 11.5378, 11.8156, 12.0933, 12.3711, 12.6489, 12.9267, 13.2044, 13.4822, 13.7600, 14.2600, 15.4533, 16.6467, 17.8400, 19.0333, 20.2267, 21.4200, 22.6133, 23.8067, 25.0000]';


lambda_i = @(lambda,beta) 1./(lambda+0.08*beta)-(0.035/(c11+c12*beta.^3));

f1 = @(lambda) c4./lambda;
f2 = @(lambda,beta) (c5*lambda_i(lambda,beta))-c6*beta-c7*beta^c8-c9;
f3 = @(lambda,beta) exp(-c10*lambda_i(lambda,beta));

dlambda_i_dlambda = @(lambda,beta) -1./(lambda+0.08*beta).^2;
dlambda_i_dbeta = @(lambda,beta) -0.08./(lambda+0.08*beta).^2+(3*0.035*c12*beta.^2)./(c11+c12*beta.^3).^2;

cQ_tilde = @(lambda,beta) c1*(1+c2*(beta+c3).^0.5)+f1(lambda).*f2(lambda,beta).*f3(lambda,beta);
cQ = @(lambda,beta) cQ_tilde(lambda,beta).*(1+sign(cQ_tilde(lambda,beta)))./2;
dcQ_dlambda = @(lambda,beta) -c4./lambda.^2.*f2(lambda,beta).*f3(lambda,beta)+f1(lambda).*(c5.*dlambda_i_dlambda(lambda,beta)).*f3(lambda,beta)+f1(lambda).*f2(lambda,beta).*(-c10.*exp(-c10*lambda_i(lambda,beta)).*dlambda_i_dlambda(lambda,beta));
dcQ_dbeta = @(lambda,beta) 0.*(c5.*dlambda_i_dbeta(lambda,beta)-c6-c7*c8*beta.^(c8-1)).*f3(lambda,beta)+f1(lambda).*f2(lambda,beta).*(-c10.*exp(-c10*lambda_i(lambda,beta)).*dlambda_i_dbeta(lambda,beta))+0.5*c1*c2*(beta+c3).^(-0.5);

lambda = @(V,omega_r) omega_r*para_mdl.R./V;

k_omega = zeros(length(v),1); % initialize the vector for k_omega
k_beta = zeros(length(v),1); % initialize the vector for k_beta
k_V = zeros(length(v),1); % initialize the vector for k_V
for i = 1:length(v)
    k_omega(i) = 0.5*para_mdl.rho*v(i)^2*pi*para_mdl.R^3*(para_mdl.R/v(i)*dcQ_dlambda(lambda(v(i),omega_r),beta(i))); % calculate k_omega for the i-th wind speed
    k_beta(i) = 0.5*para_mdl.rho*v(i)^2*pi*para_mdl.R^3*dcQ_dbeta(lambda(v(i),omega_r),beta(i)); % calculate k_beta for the i-th wind speed
    k_V(i) = 0.5*para_mdl.rho*pi*para_mdl.R^3*(2*v(i)*cQ(lambda(v(i),omega_r),beta(i))-omega_r*para_mdl.R*dcQ_dlambda(lambda(v(i),omega_r),beta(i))); % calculate k_V for the i-th wind speed
end

% % Create a table with the results
% T = table(v, k_omega, k_beta, k_V, 'VariableNames', {'WindSpeed', 'k_omega', 'k_beta', 'k_V'});

tau_ref = 3.5;
a = -k_omega/J;
b = k_beta/J;
k_p = 1./(b.*tau_ref);
k_I = a./(b.*tau_ref);


omega_g_R1_max = 0.6 * para_mdl.omega_g_R;
omega_g_R1p5_max = 0.75 * para_mdl.omega_g_R;
omega_g_R2_max = 0.9 * para_mdl.omega_g_R;
omega_g_R2p5_max = 0.6 * para_mdl.omega_g_R;
omega_g_R = para_mdl.omega_g_R;
Tg_R = para_mdl.Tg_R;



