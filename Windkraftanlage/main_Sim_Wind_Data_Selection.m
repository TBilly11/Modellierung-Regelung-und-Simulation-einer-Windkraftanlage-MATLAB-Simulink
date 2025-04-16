%% 
%% 
%-------------------------------------------------------------------------------
% main_Sim_Wind_Data_Selection                                   22.11.2020
%-------------------------------------------------------------------------------
%clear all
                                                                    

%% load wind data
cd 01_Wind

% gust_data.mat               
% gust_data_single_gusts.mat  
load gust_data.mat               

% turb_data_200s.mat          
% turb_data_highwind.mat      
% TurbWind_12ms_grid25.mat    
% turb_data.mat        
% turb_data_lowwind.mat  
load turb_data.mat                  
cd ..


%% wind data selection 

% wind type selection
para_mdl.WindInputType                  = 2;                        % 1: constant wind speed, 2: IECwind  Wind gust, 3: Turbulenter Wind,   4: wind speed steps 

para_mdl.v                              = 3; % 11.26; %18;%              % constant wind speed value if WindInputType = 1 was selected

para_mdl.v_mean                         = 18;                       % average wind speed as an input for the disturbance observer 

% select initial wind speed for IEC gust
gust_nr                                 = 7;                        %  1: 6 m/s, 2: 8 m/s, 3: 10 m/s, 4: 11 m/s, 5: 12 m/s, 6: 14 m/s, 7: 16 m/s, 8: 18 m/s
gust_wind_data.time                     = gust_data(:,1);
gust_wind_data.signals.values           = gust_data(:,gust_nr+1);
gust_wind_data.signals.dimensions       = 1;

% turbulent wind
turb_nr                                 = 4;                        % selection of the average wind speed for turbulent wind 
                                                                    % 1: 6, 2: 8, 3: 9, 4: 10, 5: 11, 6: 12, 7: 14, 8: 16, 9: 18, 10: 20, 11: 22, 12: 24                                                                   % 24
turb_wind_data.time                     = turb_data(:,1);
turb_wind_data.signals.values           = turb_data(:,turb_nr+1);
turb_wind_data.signals.dimensions       = 1;



%% open and start simulation 
dT = 0.01;   % sampling time
Tmax = 100;  % total simulation time
soptions = simset('Solver','ode4','FixedStep', dT);
open('sim_wind_speed_generator');
tic
sim('sim_wind_speed_generator', Tmax, soptions);
toc

%% plot results for documentation (book)
% ...

%% save extended para_mdl structure
%save para_mdl_sim para_mdl


