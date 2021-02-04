%% ROCKET LAUNCH

clear all;
close all;
clc;

%% 1. Definition of Constants, Parameters and Variables

% 1.1. CONSTANTS
    
    G = 6.67408e-11;            % Gravitational Constant                [m^3/(kg*s^2)]
    R = 8.31432;                % Universal Constant for Ideal Gases    [J/mole*K]
    
    % 1.1.1. Earth
    R_Earth = 6371.0e3;         % Radius of Earth                       [m]
    M_Earth = 5.9724e24;        % Earth's Mass                          [kg]
    g0_Earth = 9.80665;         % Acceleration at Earth's surface       [m/s^2]
    T0_Earth = 288.15;          % US Standard Sea Level Temperature     [K]
    P0_Earth = 101325;          % Pressure at Sea Level                 [Pa]
    Mm_Earth = 28.9644*10^-3;   % Molecular Mass                        [kg*mole^-1]
    gamma_gas_Earth = 1.4;      % Earth's air specific heats relation   [adim]
    
    % Earth's atmospheric layers altitude           [m]
    H_layer_Earth = 1e3*[0 11 20 32 47 52 61 69 79 90 100 110 117.776];    
    % Earth's atmopheric layers thermal gradient    [K/m]    
    lambda_layer_Earth = 1e-3*[-6.5 0 1 2.8 0 -2 -4 -3 0 2 4.36 16.4596 0];
    
    % 2.3. Compute base temperatures and pressures
    [Tb, Pb] = getBaseTemperaturePressure(R, g0_Earth, T0_Earth, P0_Earth, Mm_Earth, H_layer_Earth, lambda_layer_Earth);

    
% 1. 3. VARIABLES
syms V;                     % Vehicle velocity                          [m/s]
syms gamma;                 % Vehicle flight path angle                 [rad]
syms h;                     % Vehicle altitude                          [m]
syms r;                     % Vehicle flight range                      [m]
syms t;                     % Vehicle time of flight                    [s]
syms g;                     % Acceleration                              [m/s^2]
syms Q;                     % Dynamic pressure at h                     [Pa]
syms W;                     % Vehicle weight                            [N]
syms C_D;                   % Drag Coefficient                          [adim]
syms A;                     % Vehicle reference area used in C_D        [m^2]
%syms beta;                    % Ballistic coefficient                  [N/m^2]
syms C_L;                   % Veicle lift coefficient                   [adim]



%% Input data and Previous calculation
                                               
Pc = 100 * P0_Earth;                                             % Pa
Tc = 3500;                                                       % K
gamma = 1.25;
m0 = 540e3;                                                      % kg
mW = 16e-3;                                                      % kg/mol
m_dot = 2000;                                                    % kg/sec
Ae_At = (2:0.5:80)';
R_comb = R/(mW);
c_car = sqrt(R_comb*Tc)/0.6581;                                  % 0.6581 from tables 
At = (m_dot*c_car)/Pc;
MPF_e(:,1) = 0.6581./Ae_At(:,1);

% Solving MFP equation to obtain Mach Number
Me = zeros(length(MPF_e),1);
for i=1:length(MPF_e)
    syms M;
    M = vpasolve(MPF_e(i,1) == sqrt(gamma)* M / ((1 + (gamma-1)/2 *M^2)^((gamma+1)/(2*(gamma-1)))) , M, [1,6]);
    Me(i,1) = M;
end

P_e = zeros(length(Ae_At),1);
for t=1:length(Ae_At)
    P_e (t,1) = Pc/(1+Me(t,1)^2*((gamma-1)/2))^(gamma/(gamma-1));
end

%% Numerical integration

final_time = 120;                    % Final time for computations   [s]
time_step = 0.5;                    % Step size                     [s]
time = 0:time_step:final_time;      % Calculates upto final_time 

h = zeros(length(Ae_At),length(time));
v = zeros(length(Ae_At),length(time));

h(:,1) = 0;                                                      % initial condition of altitude
v(:,1) = 0;                                                      % initial condition of velocity
 
F_th = @(t,r,s) s;                                                                 % Altitude EDO

F_tv = @(t,r,s,p,q,c,m,ratio) (Pc*At*(((2/(gamma+1))^((gamma+1)/(2*(gamma-1))))*(gamma*m+(1/m))/(sqrt(1+(m^2*(gamma-1)/2)))-((p/Pc)*ratio))/ (m0-m_dot*t)) - g0_Earth;    % Velocity EDO

for i=1:length(Ae_At)
    
    for j=1:(length(time)-1)                                     % calculation loop

        [T, P, rho, a] = computeAtmosphere(Tb, Pb, H_layer_Earth, lambda_layer_Earth, R, g0_Earth, Mm_Earth, h(i,j));
        
        if P_e(i,1) < (0.4*P)
          
        Me_prima(i,1) = sqrt((2/(gamma-1)) * ((Pc/(0.4*P))^((gamma-1)/gamma)-1));
        Ae_At_prima(i,1) = 0.6581 / (sqrt(gamma)*Me_prima(i,1)/(1+((gamma-1)/2)*Me_prima(i,1)^2)^((gamma+1)/(2*(gamma-1))));
        
        k_1h = F_th(time(j),h(i,j),v(i,j));
        k_1v = F_tv(time(j),h(i,j),v(i,j),P,rho,a,Me_prima(i,1),Ae_At_prima(i,1));

        k_2h = F_th(time(j)+0.5*time_step,h(i,j)+0.5*time_step*k_1h,v(i,j)+0.5*time_step*k_1v);
        k_2v = F_tv(time(j)+0.5*time_step,h(i,j)+0.5*time_step*k_1h,v(i,j)+0.5*time_step*k_1v,P,rho,a,Me_prima(i,1),Ae_At_prima(i,1));

        k_3h = F_th((time(j)+0.5*time_step),(h(i,j)+0.5*time_step*k_2h),(v(i,j)+0.5*time_step*k_2v));
        k_3v = F_tv((time(j)+0.5*time_step),(h(i,j)+0.5*time_step*k_2h),(v(i,j)+0.5*time_step*k_2v),P,rho,a,Me_prima(i,1),Ae_At_prima(i,1));

        k_4h = F_th((time(j)+time_step),(h(i,j)+k_3h*time_step),(v(i,j)+k_3v*time_step));
        k_4v = F_tv((time(j)+time_step),(h(i,j)+k_3h*time_step),(v(i,j)+k_3v*time_step),P,rho,a,Me_prima(i,1),Ae_At_prima(i,1));

             
         %Cf_v (i,1) = ((2/(gama+1))^((gama+1)/(2*(gama-1))))*(gama*Me(i,1)+(1/Me(i,1)))/(sqrt(1+(Me(i,1)^2*(gama-1)/2)));
             
         else 
             
        k_1h = F_th(time(j),h(i,j),v(i,j));
        k_1v = F_tv(time(j),h(i,j),v(i,j),P,rho,a,Me(i,1),Ae_At(i,1));

        k_2h = F_th(time(j)+0.5*time_step,h(i,j)+0.5*time_step*k_1h,v(i,j)+0.5*time_step*k_1v);
        k_2v = F_tv(time(j)+0.5*time_step,h(i,j)+0.5*time_step*k_1h,v(i,j)+0.5*time_step*k_1v,P,rho,a,Me(i,1),Ae_At(i,1));

        k_3h = F_th((time(j)+0.5*time_step),(h(i,j)+0.5*time_step*k_2h),(v(i,j)+0.5*time_step*k_2v));
        k_3v = F_tv((time(j)+0.5*time_step),(h(i,j)+0.5*time_step*k_2h),(v(i,j)+0.5*time_step*k_2v),P,rho,a,Me(i,1),Ae_At(i,1));

        k_4h = F_th((time(j)+time_step),(h(i,j)+k_3h*time_step),(v(i,j)+k_3v*time_step));
        k_4v = F_tv((time(j)+time_step),(h(i,j)+k_3h*time_step),(v(i,j)+k_3v*time_step),P,rho,a,Me(i,1),Ae_At(i,1));   
             
         end
       
        h(i,j+1) = h(i,j) + (1/6)*(k_1h+2*k_2h+2*k_3h+k_4h)*time_step;  % main altitude equation
        v(i,j+1) = v(i,j) + (1/6)*(k_1v+2*k_2v+2*k_3v+k_4v)*time_step;  % main altitude equation

    end

end


figure(1);
hold on;
str = sprintf("\\textbf{Final altitude vs. $A_e / A_t$ at $t_b = %d \\ \\mathrm{s}$ (Without drag)}", final_time);
title(str);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(Ae_At, h(:,end)/1e3, 'b');
xlabel("$A_e / A_t$");
ylabel("Altitude $\left( \mathrm{km} \right)$");
set(gcf, 'units', 'centimeters', 'position', [1, 1, 18, 15]);
grid on;
grid minor;
box on;
hold off;

figure(2);
hold on;
str = sprintf("\\textbf{Final velocity vs. $A_e / A_t$ at $t_b = %d \\ \\mathrm{s}$ (Without drag)}", final_time);
title(str);
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(Ae_At, v(:,end), 'b');
xlabel("$A_e / A_t$");
ylabel("Velocity $\left( m / s \right)$");
set(gcf, 'units', 'centimeters', 'position', [19, 1, 18, 15]);
grid on;
grid minor;
box on;
hold off;


