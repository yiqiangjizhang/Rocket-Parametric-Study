function [Mach_menor,Mach_mayor] = MPF_function(gamma,MPF)
    %Devuelve dos valors de Mach, el valor pequeño y el grande
    %Para casos concretos de MPF pueden coincidir
    Mach_menor=vpasolve(MPF == sqrt(gamma)* M / ((1 + (gamma-1)/2 *M^2)^((gamma+1)/(2*(gamma-1)))));
    Mach_mayor=vpasolve(MPF == sqrt(gamma)* M / ((1 + (gamma-1)/2 *M^2)^((gamma+1)/(2*(gamma-1)))) , M,'Random',true);
end

