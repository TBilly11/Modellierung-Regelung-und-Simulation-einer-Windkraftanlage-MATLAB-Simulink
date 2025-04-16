%----------------------------------------------------------------------------------
% calc_CQ_CT_BEM_v6.m                                      HS   20/03/2022
%----------------------------------------------------------------------------------
% Description: Calculation of the CT and CQ look up table using the blade
%              element method. 
%
%              This requires the following data: 
%
%              1. Description of the blade elements defined by    
%                  dr: blade element length in [m}
%                  r : blade element distance to the axis of rotation [m]
%                  t : blade chord in [m]
%                  alpha_twist: angle of each balde element related to
%                      the rotor plane
%                  NFoil: ID for the allocation of an airfoil profile from blade profile data
%
%               2. Blade profile data: 
%                       - each blade element has an airfoil shape with constant cross section
%                       - airfoil profiles are mapped to the blade element via IDs 
%                       - drag and lift coefficients of the air foils must be available for  
%                         all occurring angles of attack
%              
%                3. Turbine parameter 
%                       R   rotor radius 
%                       z   number of blades 
%                       v1  design wind speed in [m/s]
%-------------------------------------------------------------------------------------

%clear all;

%-----------------------------
% 1. Load and set parameter
%-----------------------------
NREL_ON = 1; % flag for switch ON / OFF exception handling for airfoil #1 and #2
GGS_ON  = 0; 

if NREL_ON==1
    disp('Load and set parameter for NREL %MW turbine')
    load para_mdl_WT.mat             % turbine parameter: R, number of blades, air density  
    load Foils_NREL_5MW.mat          % airfoile data with cl(alpha) and cd(alpha)
    load BladeSec_NREL_5MW_Data.mat  % blade element data
    noOfBladeElements = length(BladeSec_Data); % number of  blade elements

    % Turbine parameter 
    R       = para_mdl.R;      % rotor radius [m] 
    z       = para_mdl.N;      % number of blades 
    v1      = 8;               % Design wind speed [m/s]
    rho     = para_mdl.rho;    % air density [kg/m^3] @ standard conditions

    % lambda_d = para_mdl.R * ??? / v1  % Design tip-speed ratio
else
 if GGS_ON==1
    disp('Load and set parameter for GGS turbine')
    load Foils_GGS.mat           % airfoile data with cl(alpha) and cd(alpha)
    load BladeSec_GGS_Data.mat   % blade element data
    noOfBladeElements = length(BladeSec_Data); % number of  blade elements

    % Turbine parameter 
    R       = 55 + 2.2275 ;  % rotor radius [m]
    z       = 3;          % number of blades  
    v1      = 8;          % design wind speed [m/s]
    rho     = 1.2250;     % air density [kg/m^3] @ standard conditions

 else
    disp('Choose suitable data base either NREL_ON or GGS_ON') 
 end
end
    

%-----------------------------------------
% 2. Load OR set lambda and beta array
%-----------------------------------------------------------------------

% load lambda and beta array
%load lambda_beta_array.mat      

disp('Set lambda and beta array')

% set lambda array
lambda_max = 17.5;
lambda_array = linspace(0,lambda_max,71);

% set beta array
beta_grad_max = 90;

% rough mesh
beta_array = linspace(0,beta_grad_max,19)*pi/180;  % {0, 5, ... , 85, 90}

% medium mesh
% beta_array = linspace(0,beta_grad_max,46)*pi/180;     % {0, 2, 4, ...  86, 88, 90} 

% fine mesh
% beta_array = linspace(0,beta_grad_max,91)*pi/180;     % {0, 1, 2, ...  88, 89, 90} 

% very fine mesh
%beta_array = linspace(0,beta_grad_max,181)*pi/180;     % {0, 0.5, 1, ...  89, 89.5, 90} 


%-----------------------------------------------------------------------
% 3. Calc CQ and CT with BEM algorithm according to Gasch/Twele 
%-----------------------------------------------------------------------
disp('Calculation of CQ and CT with BEM algorithm')

omr = lambda_array * v1/R; %  % rotor speed array for all lambda   

BladeSec_Data(end,2) = 0.15;   % set twist angle at the tip (last blade element) to 0.15
                               % caused an increase of CPmax = 0.4615 to 0.4661 

BlSec.r           = BladeSec_Data(:,1);   % blade element distance to rotation axis
BlSec.alpha_twist = BladeSec_Data(:,2);   % aero twist [grad]
BlSec.dr          = BladeSec_Data(:,3);   % blade element length
BlSec.t           = BladeSec_Data(:,4);   % blade chord
BlSec.NFoil       = BladeSec_Data(:,5);   % airfoil ID (number) 


for ibeta = 1:length(beta_array)    % loop over all beta values    
  
  ilambda = 0;
  beta_array(ibeta)*180/pi  % beta in degrees
  
  for l = 1:length(lambda_array)       % loop over all lambda values

        ilambda = ilambda + 1;
        om = omr(l);
        beta = beta_array(ibeta);  % [rad]

        for i = 1:noOfBladeElements       % loop over all blade elements 

            r           = BlSec.r(i);
            dr          = BlSec.dr(i);
            t           = BlSec.t(i);
            alpha_twist = BlSec.alpha_twist(i) * pi/180 + beta;   % [grad] -> [rad], pitch angle as twist offset 
            NFoil       = BlSec.NFoil(i);

            f_old = 0;
            alpha_decr = 1 * pi/180;


            % BEM algorithm according to Gasch/Twele, Section 6.9 

            fmin = 0.000001;      % stop criterion of all iterations 
            NrOfIterations_Max = 1000; % stop criterion of all iterations

            % 1. Initial value: alpha = alpha1

            c1 = sqrt(v1^2 + (om * r)^2);    % magnitude of inflow speed in front of the rotor (undisturbed)

            alpha1 = acos(om*r/c1);          % angle between rotational speed and inflow speed [rad]
            alpha  = alpha1;  


            % 2. Limits of the angle of attack: alpha_min < alpha < alpha_max 

            sin_alpha_max = z * sqrt(1 - (r/R)^2) / (2*pi * r/R);
            sin_alpha_min = sin(2/3*alpha1);

            stop = false;
            NrOfIterations(i) = 0;

            
            % Iteration loop for calculation of alpha where alphaA = alpha - twist
            while stop == false

                    NrOfIterations(i) = NrOfIterations(i) + 1;

                    % 3. Calculation of the lift and drag coefficient

                    alpha_A = alpha - alpha_twist;     % angle of attack [rad]

                    if ((NFoil <= 2)  &&  (NREL_ON==1))
                        Cl = Foil(NFoil).Cl;        % airfoils 1 and 2 are cylindrical, contains only one value for Cl and Cd
                        Cd = Foil(NFoil).Cd;
                    else
                        Cl = interp1(Foil(NFoil).alpha, Foil(NFoil).Cl, alpha_A * 180/pi, 'linear', 0);  % Cl(alpha_A * 180/pi)
                        Cd = interp1(Foil(NFoil).alpha, Foil(NFoil).Cd, alpha_A * 180/pi, 'linear', 0);  % Cd(alpha_A * 180/pi)
                
                    end

                    sin_alpha = sin(alpha);

                    x = sin_alpha;


                    % 4. Corrections according to Glauert OR Prandtl  

                    if      x < sin_alpha_min

                        y = sin_alpha/(sin_alpha_min);
                        x = 1/4 * sin_alpha_min * sqrt(9 - 2*y^2 + 9*y^4);      % Corrections according to Glauert

                    elseif  x > sin_alpha_max

                        x = sin_alpha_max;                                      % Corrections according to Prandtl

                    end


                    % 5. Iterative equation for calculation of alpha which fulfill f(alpha) = 0 

                    f(i) = Cl - (8*pi*r/(z*t)*x + Cd) * tan(alpha1 - alpha);


                    % 6. Check and terminate condition

                    if      f(i) < 0

                        if f_old > 0
                            alpha_decr = alpha_decr - 0.5 *alpha_decr;
                        end

                        alpha = alpha + alpha_decr;

                    elseif  f(i) > 0

                        if f_old < 0
                            alpha_decr = alpha_decr - 0.5 *alpha_decr;
                        end

                        alpha = alpha - alpha_decr;
                    end

                    % terminate condition
                    if ((abs(f(i)) < fmin)  || (NrOfIterations(i) > NrOfIterations_Max))  
                        stop = true; 
                    end

                    f_old = f(i);

            end     % while (Iteration)  -> TRUE alpha


            % 7. Calculation of the inflow velocity within the rotor plane

            c = c1 * cos(alpha1 - alpha) * (8*pi*r/z * x)/(8*pi*r/z * x + t * Cd);


            % 8. Forces at blade element

            dA(i) = rho/2 * c^2 * t * dr * Cl;
            dW(i) = rho/2 * c^2 * t * dr * Cd;

            dS(i) = dA(i) * cos(alpha) + dW(i) * sin(alpha);
            dU(i) = dA(i) * sin(alpha) - dW(i) * cos(alpha);

            dM(i) = dU(i) * r;

         end   % loop over all blade elements
         
         % Calculation of trust force S and rotor torque M as sum of all blade elements  
         S(ilambda, ibeta) = 0;   % init values
         U(ilambda, ibeta) = 0;
         M(ilambda, ibeta) = 0;
         
         for i = 1:noOfBladeElements
        S(ilambda, ibeta) = S(ilambda, ibeta) + z * dS(i);
        U(ilambda, ibeta) = U(ilambda, ibeta) + z * dU(i);
        M(ilambda, ibeta) = M(ilambda, ibeta) + z * dM(i);
         end
         
         % Calculation of CT and CQ using FT and TR of the current tube theory
         CT(ilambda, ibeta) = S(ilambda, ibeta) / (0.5*rho*pi*R^2*v1^2);
         CQ(ilambda, ibeta) = M(ilambda, ibeta) / (0.5*rho*pi*R^3*v1^2);

         if CT(ilambda, ibeta) < 0
             CT(ilambda, ibeta) = 0;
         end
         if CQ(ilambda, ibeta) < 0
             CQ(ilambda, ibeta) = 0;
         end
 end  % loop over all lambda values
end  % loop over all beta values


%-----------------------------------------------------------------------
% 4. Calc CP
%-----------------------------------------------------------------------
% Calculation of the power coefficients CP= CQ * lambda

disp('Calculation of power coefficients CP = CQ * lambda')
for ibeta = 1:length(beta_array)
     CP(:,ibeta) = CQ(:,ibeta)' .* lambda_array;
end


%-----------------------------------------------------------------------
% 5. Save all coefficients  
%     CT Thrust coefficients
%     CQ Torque coefficients
%     CP Power coefficients
%-----------------------------------------------------------------------
disp('Save all coefficients')

% copy all in a structure and save in mat
Cx_maps.CT = CT;
Cx_maps.CQ = CQ;
Cx_maps.CP = CP;
Cx_maps.lambda_array = lambda_array;
Cx_maps.beta_array   = beta_array;


%save Cx_maps_very_fine_v5.mat CT CQ CP lambda_array beta_array
%save Cx_maps_very_fine_v5.mat Cx_maps
save Cx_maps.mat Cx_maps




