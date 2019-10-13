%calculate the overall stress
function [sigma_collagen_lp, sigma_collagen_dsm, sigma_elastin, sigma] = stress_cal(DR,LR,K,KE,lam)
k_elastin = KE;
k_collagen = K;
%k_muscle_p = K(4);
ratio = 4;
% Lamina propria (LP) collagen attachment stretch distribution
col_rec_max_lp = LR(2); %maximum attachment stretch 1
col_rec_min_lp = LR(1); %minimum attachment stretch
col_rec_mod_lp = LR(3); %mode attachment stretch
col_rec_lp_width = col_rec_max_lp-col_rec_min_lp; %distribution width 0.2
col_rec_lp_skew = (col_rec_mod_lp-col_rec_min_lp)/col_rec_lp_width; %distribution skew 0.4
% Destrusor (DSM) collagen attachment stretch
% Lamina propria (LP) collagen attachment stretch distribution
col_rec_max_dsm = DR(2); %maximum attachment stretch 1
col_rec_min_dsm = DR(1); %minimum attachment stretch
col_rec_mod_dsm = DR(3); %mode attachment stretch
col_rec_dsm_width = col_rec_max_dsm-col_rec_min_dsm; %distribution width 0.2
col_rec_dsm_skew = (col_rec_mod_dsm-col_rec_min_dsm)/col_rec_dsm_width; %distribution skew 0.4
% % Smooth muscle attachment stretch
% smc_att_dsm = 1.52;
% % Smooth muscle recruitment stretch at initial state
%  smc_rec_dsm = 1.3;
% %% Define the muscle stress function (passive)
% % Define the smooth muscle stretch (= attachment stretch at initial state)
% lambda_smc  = lam/smc_rec_dsm;
% % Muscle stress function (passive)
% sigma_muscle_p = lambda_smc.^2 .* k_muscle_p .* (1-(1./(lambda_smc.^6)));
% %% Define the elastin stress function
  sigma_elastin  = lam.^2 .* k_elastin .* (1 - (1./lam.^6));
%% Define the collagen stress function (detrusor)
% The collagen distribution (min,mod,max)
v_a_dsm = col_rec_min_dsm;
v_c_dsm = col_rec_mod_dsm;
v_b_dsm = col_rec_max_dsm;
% Gamma and delta are two common factors 
v_gamma_dsm   = k_collagen / ((v_b_dsm - v_a_dsm) * (v_c_dsm - v_a_dsm));
v_delta_dsm   = k_collagen / ((v_b_dsm - v_a_dsm) * (v_b_dsm - v_c_dsm));
x = lam;
% The collagen stress distribution (triangular PDF)            
sigma_collagen_dsm_0      =  x * 0;
sigma_collagen_dsm_ac     =  x .* v_gamma_dsm .* 2 .* ( (x + v_a_dsm) .* log(x./v_a_dsm) + 2.*(v_a_dsm - x) ) ;
sigma_collagen_dsm_cb     =  x .* v_gamma_dsm .* 2 .* ( (x + v_a_dsm).*log(v_c_dsm./v_a_dsm) + v_a_dsm - v_c_dsm + ((v_a_dsm - v_c_dsm) ./ v_c_dsm) .* x) ...
                - x .* v_delta_dsm .* 2 .* ((x + v_b_dsm).*log(x./v_c_dsm) + v_b_dsm + v_c_dsm - ((v_b_dsm + v_c_dsm) ./ v_c_dsm) .* x );
sigma_collagen_dsm_b      =  x .* v_gamma_dsm .* 2 .* ( (x + v_a_dsm).*log(v_c_dsm./v_a_dsm) + v_a_dsm - v_c_dsm + ((v_a_dsm - v_c_dsm) ./ v_c_dsm) .* x) ...
                - x .* v_delta_dsm .* 2 .* ((x + v_b_dsm).*log(v_b_dsm./v_c_dsm) - v_b_dsm + v_c_dsm - ((v_b_dsm - v_c_dsm) ./ v_c_dsm) .* x);
% The function of collagen stress          
sigma_collagen_dsm       = sigma_collagen_dsm_0.*(x<v_a_dsm)...
                + sigma_collagen_dsm_ac.*( x>=v_a_dsm & x<v_c_dsm)...
                + sigma_collagen_dsm_cb.*(x>=v_c_dsm & x<=v_b_dsm)...
                + sigma_collagen_dsm_b.*(x>v_b_dsm);
            %% Define the collagen stress function (lamina propria)
% The collagen distribution (min,mod,max)
v_a_lp = col_rec_min_lp;
v_c_lp = col_rec_mod_lp;
v_b_lp = col_rec_max_lp;
% Gamma and delta are two common factors 
v_gamma_lp   =   ratio*k_collagen / ((v_b_lp - v_a_lp) * (v_c_lp - v_a_lp));
v_delta_lp   =   ratio*k_collagen / ((v_b_lp - v_a_lp) * (v_b_lp - v_c_lp));
% The collagen stress distribution (triangular PDF)            
sigma_collagen_lp_0      = x * 0;
sigma_collagen_lp_ac     = x .* v_gamma_lp .* 2 .* ( (x + v_a_lp) .* log(x./v_a_lp) + 2.*(v_a_lp - x) ) ;
sigma_collagen_lp_cb     = x .* v_gamma_lp .* 2 .* ( (x + v_a_lp).*log(v_c_lp./v_a_lp) + v_a_lp - v_c_lp + ((v_a_lp - v_c_lp) ./ v_c_lp) .* x) ...
                - x .* v_delta_lp .* 2 .* ((x + v_b_lp).*log(x./v_c_lp) + v_b_lp + v_c_lp - ((v_b_lp + v_c_lp) ./ v_c_lp) .* x );
sigma_collagen_lp_b      =  x .* v_gamma_lp .* 2 .* ( (x + v_a_lp).*log(v_c_lp./v_a_lp) + v_a_lp - v_c_lp + ((v_a_lp - v_c_lp) ./ v_c_lp) .* x) ...
                - x .* v_delta_lp .* 2 .* ((x + v_b_lp).*log(v_b_lp./v_c_lp) - v_b_lp + v_c_lp - ((v_b_lp - v_c_lp) ./ v_c_lp) .* x);
% The function of collagen stress          
sigma_collagen_lp       = sigma_collagen_lp_0.*(x<v_a_lp)...
                + sigma_collagen_lp_ac.*( x>=v_a_lp & x<v_c_lp)...
                + sigma_collagen_lp_cb.*(x>=v_c_lp & x<=v_b_lp)...
                + sigma_collagen_lp_b.*(x>v_b_lp);
%% Define the total stress
sigma = sigma_collagen_lp+sigma_collagen_dsm+sigma_elastin; 
% +sigma_elastin+sigma_muscle_p;
end