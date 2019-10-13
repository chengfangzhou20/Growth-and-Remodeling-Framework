% Collagen remodeling: update the recruitment distribution
function [UR,UM] = collagen_remodeling(M,R,At,Alpha,Beta,lambda)
global deltaT
% UR: updated recruitment distribution; R: current recruitment
% distribution; At: attachment distribution;
% Calculate collagen distribution
[C,W,S] = collagen_stretch(R,lambda);
WR = R(2)-R(1);
SR = (R(3)-R(1))/WR;
% Collagen recruitment stretch remoldeling 
UR(2) = R(2)+ Alpha*(C(1)- At(1))/(At(1))*deltaT;
UR(1) = UR(2)-WR;
UR(3) = UR(1)+SR*WR;
% UR(1) = R(1)+ Alpha*(C(2)- At(2))/(At(2))*deltaT;
% UR(3) = R(3)+ Alpha*(C(3)- At(3))/(At(3))*deltaT;
% Collagen growth (mass density): controlled by the recruitment maximum stretch
UM = M + Beta*(C(2)- At(2))/(At(2))*deltaT;
end