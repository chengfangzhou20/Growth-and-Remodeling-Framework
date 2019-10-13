
%% Define the collagen stress function
function sigma_collagen = collagen_stress(x,R,ratio) %tissue stretch x
%Att is the attachment stretch, Att(1) min, Att(2) max, Att(3) mod,
global k_collagen % Global variable  k_collagen(collagen stiffness)
% The collagen distribution (min,mod,max)
v_a = R(1); % min recruitment stretch
v_c = R(3); % mod recruitment stretch
v_b = R(2); % max recruitment stretch
% Gamma and delta are two common factors 
v_gamma   = ratio*k_collagen / ((v_b - v_a) * (v_c - v_a));
v_delta   = ratio*k_collagen / ((v_b - v_a) * (v_b - v_c));
% The collagen stress distribution (triangular PDF)            
sigma_collagen_0      =  x * 0;
sigma_collagen_ac     =  x * v_gamma * 2 * ( (x + v_a) * log(x/v_a) + 2*(v_a - x) ) ;
sigma_collagen_cb     =  x * v_gamma * 2 * ( (x + v_a)*log(v_c/v_a) + v_a - v_c + ((v_a - v_c) / v_c) * x) ...
                - x * v_delta * 2 * ((x + v_b)*log(x/v_c) + v_b + v_c - ((v_b + v_c) / v_c) * x );
sigma_collagen_b      =  x * v_gamma * 2 * ( (x + v_a)*log(v_c/v_a) + v_a - v_c + ((v_a - v_c) / v_c) * x) ...
                - x * v_delta * 2 * ((x + v_b)*log(v_b/v_c) - v_b + v_c - ((v_b - v_c) / v_c) * x);
% The function of collagen stress          
sigma_collagen       = sigma_collagen_0.*(x<v_a)...
                + sigma_collagen_ac.*( x>=v_a & x<v_c)...
                + sigma_collagen_cb.*(x>=v_c & x<=v_b)...
                + sigma_collagen_b.*(x>v_b);
end