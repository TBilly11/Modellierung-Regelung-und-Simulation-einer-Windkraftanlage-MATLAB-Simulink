% Define state enumeration
States = {'Region_1', 'Region_1_5', 'Region_2', 'Region_2_5', 'Region_3'};

% Define parameters
wg_R1_max = 10; % replace with actual value
wg_R1_5_max = 20; % replace with actual value
wg_rated = 30; % replace with actual value
dwg = 5; % replace with actual value
Tg_max = 100; % replace with actual value

% Initialize state
current_state = States{1}; % Start with the first state

% Simulate transitions (example values, replace with actual inputs)
Tg = 90; % replace with actual value
wg = 15; % replace with actual value
b = 1; % replace with actual value

% Simulation time steps
for time_step = 1:10
    % Evaluate transition conditions
    switch current_state
        case 'Region_1'
            if wg > wg_R1_max
                current_state = States{2};
            end
            
        case 'Region_1_5'
            if wg > wg_R1_5_max
                current_state = States{3};
            end
            
        case 'Region_2'
            if wg > (wg_rated - dwg)
                current_state = States{4};
            end
            
        case 'Region_2_5'
            if wg > wg_rated && Tg > Tg_max
                current_state = States{5};
            end
            
        case 'Region_3'
            % Define transition condition
            % e) Transition Region_3 ----> Region_2.5
            if wg < wg_rated && Tg < Tg_max % Replace with actual condition
                current_state = States{4};
            end
            
        case 'Region_2_5'
            % Define transition condition
            % f) Transition Region_2.5 ----> Region_2
            if wg < (wg_rated - dwg) % Replace with actual condition
                current_state = States{3};
            end
            
        case 'Region_2'
            % Define transition condition
            % g) Transition Region_2 ----> Region_1.5
            if wg < wg_R1_5_max  % Replace with actual condition
                current_state = States{2};
            end
            
        case 'Region_1_5'
            % Define transition condition
            % h) Transition Region_1.5 ----> Region_1
            if wg < wg_R1_max % Replace with actual condition
                current_state = States{1};
            end
    end
    
    % Display current state (for simulation purposes)
    disp(['Time Step: ', num2str(time_step), ', Current State: ', current_state]);
    
    % Update inputs (replace with actual input updates)
    Tg = Tg + 5;
    wg = wg + 2;
    b = mod(b + 1, 2);
end
