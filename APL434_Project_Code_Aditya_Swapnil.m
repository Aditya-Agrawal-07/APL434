% ==================================================
% ========= Author :- Swapnil Kashyap, Aditya Agrawal
% ========= Date and Time :- 15-04-2024, 15:53 
% ==================================================

% ==== CODE STARTS ====

clear all; close all;
%%
% ==== IMPORTANT PARAMETERS ====
% global alpha; alpha = deg2rad(40);
% global eps_l; eps_l = 0.04;
% global d; d = 0.00051;
% global l0; l0 = 0.083;
% global n; n = 12;
% global x0; x0 = 0.085;
% global K_x; K_x = 584;
% global rho; rho = 6450;
% global c_p; c_p = 836.8;
% global h_T; h_T = 70;
% global del_H; del_H = 24000;
% global T_inf; T_inf = 25;
% global theta_T; theta_T = 0.55*1e6;
% global V_in; V_in = 7;
% global E_A; E_A = 75*1e9;
% global E_M; E_M = 28*1e9;
% global C_A; C_A = 10*1e6;
% global C_M; C_M = 10*1e6;
% global A_s; A_s = 82.9;
% global A_f; A_f = 87.5;
% global M_s; M_s = 57.9;
% global M_f; M_f = 47.5;
% global r_A; r_A = 100*1e-8;
% global r_M; r_M = 80*1e-8;

alpha = deg2rad(40);
eps_l = 0.04;
d = 0.00051;
l0 = 0.083;
n = 12;
x0 = 0.085;
K_x = 584;
rho = 6450;
c_p = 836.8;
h_T = 70;
del_H = 24000;
T_inf = 25;
theta_T = 0.55*1e6;
E_A = 75*1e9;
E_M = 28*1e9;
C_A = 10*1e6;
C_M = 10*1e6;
A_s = 82.9;
A_f = 87.5;
M_s = 57.9;
M_f = 47.5;
r_A = 100*1e-8;
r_M = 80*1e-8;

% Time Params
ts = 0.0; t_heat = 8.0; t_cool = 15.0; 
steps = 1001; dt = t_cool/(steps-1); 

% Input Voltage
V_in = zeros(steps,1);
V_in(1:floor(steps*t_heat/t_cool)) = 7.0;

% ==== Setting up important parameters ====
A_cross = 0.25 * pi * d^2;
mass = rho * A_cross * l0;
A_c = pi * d * l0;
% a_A = pi / (A_f - A_s); b_A = -a_A / C_A;
% a_M = pi / (M_s - M_f); b_M = -a_M / C_M;

% ==== Setting up initial parameters ====
zhi0 = 1; zhi0_dot = 0;
sigma0 = (K_x * x0) / (n*A_cross*cos(alpha)); 
eps0 = 0;
T0 = T_inf; 
dx0 = 0; 

% ==== Setting up the derivatives
R_ohm0 = get_resistance(l0, A_cross, zhi0, r_M, r_A);
T0_dot = (V_in(1)/n)^2/(R_ohm0 * mass * c_p);
sigma0_dot = 2*1e2;
E0 = young_mod(zhi0, E_M, E_A);
eps0_dot = (sigma0_dot - theta_T * T0_dot)/E0;
dx0_dot = l0 * eps0_dot * sec(alpha);


% Setting up containers
TIME = linspace(ts, t_cool, steps);
STRESS = zeros(1, steps); STRESS(1) = sigma0;
STRESS_DOT = zeros(1, steps); STRESS_DOT(1) = sigma0_dot;
STRAIN = zeros(1, steps); STRAIN(1) = eps0;
STRAIN_DOT = zeros(1, steps); STRAIN_DOT(1) = eps0_dot;
STROKE = zeros(1, steps); STROKE(1) = dx0;
STROKE_DOT = zeros(1, steps); STROKE_DOT(1) = dx0_dot;
TEMP = zeros(1, steps); TEMP(1) = T0;
TEMP_DOT = zeros(1, steps); TEMP_DOT(1) = T0_dot;
MFRAC = zeros(1, steps); MFRAC(1) = zhi0;
MFRAC_DOT = zeros(1, steps); MFRAC_DOT(1) = zhi0_dot;
F_M = zeros(1, steps);

% ==== Forward Euler Time Stepping ====

ARG = zeros(1, steps);

% ==== Reverse transformation ====
for i = 2:steps

    % ==== Heat Tranfer Model ====
    [T, T_dot] = heat_transfer_model(TEMP(i-1), mass, c_p, V_in(i)/n, l0, A_cross, ...
                                     MFRAC(i-1), r_M, r_A, A_c, h_T, T_inf, ...
                                     del_H, MFRAC_DOT(i-1), dt, 1);
    TEMP(i) = T; TEMP_DOT(i) = T_dot;

    % ==== Phase Transformation Model ====
    [zhi, zhi_dot, arg] = phase_transformation_model(MFRAC(i-1), dt, ...
                                                1, A_f, A_s, C_A, ...
                                                0, M_f, M_s, C_M, ...
                                                TEMP(i), TEMP_DOT(i), ...
                                                STRESS(i-1), STRESS_DOT(i-1), ...
                                                0);
    MFRAC(i) = zhi; MFRAC_DOT(i) = zhi_dot;
    ARG(i) = arg;

    % === Dynamics Model ===
    [dx, dx_dot, F_m] = dynamics_model(STROKE(i-1), l0, STRAIN(i-1), n, STRESS(i-1), A_cross, ...
                                    alpha, K_x, x0, E_M, E_A, MFRAC(i), dt);
    STROKE(i) = dx; STROKE_DOT(i) = dx_dot; F_M(i) = F_m;

    % === Kinematics Model ===
    [eps, eps_dot] = kinematics_model(l0, STROKE(i), alpha, STROKE_DOT(i));
    STRAIN(i) = eps; STRAIN_DOT(i) = eps_dot;

    % === Constitutive Model ===
    [sigma, sigma_dot] = constitutive_model(STRESS(i-1), MFRAC(i), E_M, E_A, ...
                                            STRAIN_DOT(i), eps_l, MFRAC_DOT(i), ...
                                            theta_T, TEMP_DOT(i), dt);
    STRESS(i) = sigma; STRESS_DOT(i) = sigma_dot;

end

% Plotting relevant state variables against time
% figure('Position', [10, 10, 1920, 1080]);

% Plot TEMP vs TIME
% subplot(3, 2, 1);
figure('Position', [10, 10, 1920, 1080]);
plot(TIME, TEMP, 'b');
xlabel('Time');
ylabel('Temperature');
title('Temperature vs Time');

% Plot STRESS vs TIME
% subplot(3, 2, 2);
figure('Position', [10, 10, 1920, 1080]);
plot(TIME, STRESS, 'r');
xlabel('Time');
ylabel('Stress');
title('Stress vs Time');

% Plot STRAIN vs TIME
% subplot(3, 2, 3);
figure('Position', [10, 10, 1920, 1080]);
plot(TIME, STRAIN, 'g');
xlabel('Time');
ylabel('Strain');
title('Strain vs Time');

% Plot M_FRAC vs TIME
% subplot(3, 2, 4);
figure('Position', [10, 10, 1920, 1080]);
plot(TIME, MFRAC, 'm');
xlabel('Time');
ylabel('Martensite Fraction');
title('Martensite Fraction vs Time');

% Plot STROKE vs TIME
% subplot(3, 2, 5);
figure('Position', [10, 10, 1920, 1080]);
plot(TIME, abs(STROKE), 'k');
xlabel('Time');
ylabel('Stroke');
title('Stroke vs Time');

% Plot F_M VS TIME
figure('Position', [10, 10, 1920, 1080]);
plot(TIME, F_M);
xlabel('Time');
ylabel('Muscle Force');
title('Muscle Force vs Time');

% Adjust subplot spacing
% sgtitle('Sample Data Plots'); % Overall title for the figure

% find(MFRAC_DOT ~= 0)


function [zhi, zhi_dot, arg] = phase_transformation_model(zhi0, dt, ...
                                                     zhiM, A_f, A_s, C_A, ...
                                                     zhiA, M_f, M_s, C_M, ...
                                                     T0, T0_dot, ...
                                                     sigma0, sigma0_dot, ...
                                                     bool)
    
    % === Martensite to Austenite ===
    if bool == 0
        a_A = pi / (A_f - A_s); b_A = -a_A / C_A;
        % A_s_dash = A_s + s
        arg = a_A * (T0 - A_s) + b_A * sigma0;
        if arg <= 0
            zhi = zhiM; zhi_dot = 0;
            % disp('lessgo')
        elseif arg <= pi
            zhi = (zhiM / 2) * (1 + cos(arg));
            zhi_dot = (zhi - zhi0) / dt;
            % eta_sigma = (zhiM / 2) * (a_A / C_A) * sin(arg);
            % eta_T = (-zhiM / 2) * a_A * sin(arg);
            % zhi_dot = eta_sigma * sigma0_dot + eta_T * T0_dot;
            % disp('lessgo')
        else
            zhi = 0; zhi_dot = 0;
        end

    % === Austenite to Martensite ===
    else
        a_M = pi / (M_s - M_f); b_M = -a_M / C_M;        
        arg = a_M * (T0 - M_f) + b_M * sigma0;
        if arg >= pi
            zhi = zhiA; zhi_dot = 0;
        elseif arg >= 0
            zhi = ((1-zhiA) / 2) * cos(arg) + (1+zhiA) / 2;
            eta_sigma = ((1-zhiA) / 2) * (a_M / C_M) * sin(arg);
            eta_T = (-(1-zhiA) / 2) * a_M * sin(arg);
            zhi_dot = eta_sigma * sigma0_dot + eta_T * T0_dot;
        else
            zhi = 1; zhi_dot = 0;
        end
    end
    
end

function [T, T_dot] = heat_transfer_model(T0, mass, c_p, V_in, l0, A_cross, ...
                                         zhi0, r_M, r_A, A_c, h_T, T_inf, ...
                                         del_H, zhi0_dot, dt, ok)
    R_ohm = get_resistance(l0, A_cross, zhi0, r_M, r_A);
    % [V_in ^ 2 / R_ohm, A_c * h_T * (T0 - T_inf), ok * mass * del_H * zhi0_dot]
    T_dot = V_in ^ 2 / R_ohm - A_c * h_T * (T0 - T_inf) + ok * mass * del_H * zhi0_dot; % ok = 1 for M to A, ok = -1 for A to M
    T_dot = T_dot / (mass * c_p);
    T = T0 + dt * T_dot;
end

function [sigma, sigma_dot] = constitutive_model(sigma0, zhi0, E_M, E_A, ...
                                                eps0_dot, eps_l, zhi0_dot, ...
                                                theta_T, T0_dot, dt)
    E = young_mod(zhi0, E_M, E_A);
    sigma_dot = E * (eps0_dot - eps_l * zhi0_dot) + theta_T * T0_dot;
    % sigma_dot = E * (eps0_dot - eps_l * 0) + theta_T * T0_dot;
    sigma = sigma0 + dt * sigma_dot;
end

function [eps, eps_dot] = kinematics_model(l0, dx0, alpha, dx0_dot)

    l = len(l0, dx0, alpha);
    eps = 1 - l/l0;
    eps_dot = (dx0_dot * cos(alpha))/l;
end

function [dx, dx_dot, F_m] = dynamics_model(dx0, l0, eps0, n, sigma0, A_cross, ...
                                       alpha, K_x, x0, E_M, E_A, zhi0, dt)

    E = young_mod(zhi0, E_M, E_A);
    dx = l0 * (1 - eps0) * (n * sigma0 * A_cross * cos(alpha) - K_x * x0);
    dx = dx/(n * A_cross * (sigma0 * sin(alpha)^2 - E * (1 - eps0) * cos(alpha)^2));
    dx_dot = (dx - dx0)/dt;
    F_f = sigma0 * A_cross;
    F_m = muscle_force(n, F_f, alpha, K_x, x0);
end

function R_ohm = get_resistance(l, A_cross, zhi, r_M, r_A)
    R_ohm = (l / A_cross) * (zhi * r_M + (1-zhi) * r_A);
end

function output = len(l0, dx, alpha)
    output = sqrt(l0^2 + dx^2 - 2*l0*dx*cos(alpha));
end

function output = young_mod(zhi, E_M, E_A)
    output = zhi * E_M + (1 - zhi) * E_A;
end

function output = muscle_force(n, F_f, alpha, K_x, x0)
    output = n * F_f * cos(alpha) - K_x * x0;
end

% ==== CODE ENDS ====

% ==================================================
% ========= Author :- Swapnil Kashyap, Aditya Agrawal
% ========= Date and Time :- 15-04-2024, 15:53 
% ==================================================


