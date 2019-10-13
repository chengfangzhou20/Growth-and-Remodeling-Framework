%Calculate the collagen stretch distribution based on recruitment
%distribution or 'calculate the recruitment stretch based on the collagen
%stretch'
function [lambda_C,W,S] = collagen_stretch(R,lambda)
% R is the recruitment stretch: R(1)min, R(2)max, R(3)mod
% W and S are the width and skew of collagen distribution
% lambda is the tissue stretch, lambda_C is the collagen stretch
% Lamina propria (LP) collagen attachment stretch distribution
lambda_C(2) = lambda/R(1); %maximum attachment stretch 
lambda_C(1) = lambda/R(2); %minimum attachment stretch
lambda_C(3) = lambda/R(3); %mode attachment stretch
W = lambda_C(2)-lambda_C(1); %distribution width 
S = (lambda_C(3)-lambda_C(1))/W; %distribution skew

end