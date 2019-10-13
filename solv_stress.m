%Solve for tissue stretch with a given stress
function lambda_e = solv_stress(stress,R_LP,R_DSM)
global radius_tzero thickness_tzero ratio
% Define dimension term
D = @(x) 2.*thickness_tzero ./ (radius_tzero .* x.^3);
% Define collagen stress and pressure (LP)
S_C_LP = @(x) collagen_stress(x,R_LP,ratio);
P_C_LP = @(x) S_C_LP(x)*D(x);
% Define collagen stress and pressure (DSM)
S_C_DSM = @(x) collagen_stress(x,R_DSM,1);
P_C_DSM = @(x) S_C_DSM(x)*D(x);
% Define elastin stress and pressure
S_E = @(x) elastin_stress(x);
P_E = @(x) S_E(x)*D(x);
% Define the overall pressure 
P = @(x) P_C_LP(x)+P_C_DSM(x)+P_E(x);
% Define the overall stress
S = @(x) S_C_LP(x)+S_C_DSM(x)+S_E(x);
% Define the function to solve at the state of threshold pressure
fun_threshold = @(x) S(x) - stress;
%% Solve the function and calculate the stretch at threshold pressure. (the elastin stretch is the same as tissue stretch)
lambda_e = fzero(fun_threshold,1);
end