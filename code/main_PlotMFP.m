clear;
close all;
clc;

gamma = 1.1:0.1:1.4;

M = 0:0.01:10;

MFP = zeros(length(gamma), length(M));
for i = 1:length(gamma)
    MFP(i,:) = sqrt(gamma(i))*M./(1+(gamma(i)-1)/2*M.^2).^((gamma(i)+1)/(2*(gamma(i)-1)));
end

legend_str = cell(length(gamma),1);
for i = 1:length(gamma)
    legend_str(i,1) = {sprintf('$\\gamma = %.2f$', gamma(i))};
end

figure(1);
hold on;
title("\textbf{Mass flow parameter vs. Mach number}");
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
for i = 1:length(gamma)
    plot(M, MFP(i,:));
end
xlabel("Mach number");
ylabel("Mass flow parameter");
grid on;
grid minor;
box on;
set(gcf, 'units', 'centimeters', 'position', [1,1,18,15]);
legend(legend_str);
hold off;

