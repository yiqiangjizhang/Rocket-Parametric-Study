function [T, P, rho, c] = computeAtmosphere(Tb, Pb, H_layer, lambda, R, g0, Mm, H)
%--------------------------------------------------------------------------
% Inputs:
%   - Tb        Vector of base temperatures, 1 x length(H_layer)    [K]
%   - Pb        Vector of base pressure, 1 x length(H_layer)        [Pa]
%   - H_layer   Atmospheric layers altitude                         [m]
%   - lambda    Atmospheric layers thermal gradients                [K/m]
%   - gamma     Ratio of specific heats                             [adim]
%   - R         Universal gas constant                              [J/(kg K)]
%   - g0        Acceleration at planet's surface                    [m/s^2]
%   - Mm        Planet's air molar mass                             [kg/mol]
%   - H         Altitude                                            [m]
%--------------------------------------------------------------------------
% Outputs:
%   - T         Temperature at altitude H                           [K]
%   - P         Pressure at altitude H                              [Pa]
%   - rho       Density at altitude H                               [kg/m^3]
%   - c         Speed of sound at altitude H                        [m/s]
%--------------------------------------------------------------------------

% Initialize variables
T = 0;
P = 0;
rho = 0;
c = 0;

% Find atmospheric layer
found = 0;
layer = 1;
while (layer <= length(H_layer)-1) && (found == 0)
    if (H_layer(layer) <= H) && (H < H_layer(layer+1))
        found = 1;        
    else        
        layer = layer + 1;
    end
end

% Compute properties
if found == 1    
    % Temperature
    T = Tb(layer) + lambda(layer)*(H - H_layer(layer)); 
    % Pressure
    if lambda(layer) == 0
        P = Pb(layer)*exp(-g0*Mm*(H-H_layer(layer))/(R*Tb(layer)));
    else
        P = Pb(layer)*((Tb(layer)/(Tb(layer) + lambda(layer)*(H-H_layer(layer))))^(g0*Mm/(R*lambda(layer))));
    end
    % Density
    rho = P/((R/Mm)*T);
    % Specific heats and ratio of specific heats at constant pressure
    cp = 1034.09 - 2.849e-1*T + 7.817e-4*T^2 - 4.971e-7*T^3 + 1.077e-10*T^4;
    cv = cp - R/Mm;
    gamma = cp/cv;
    % Speed of sound
    c = sqrt(gamma*(R/Mm)*T);
end

end