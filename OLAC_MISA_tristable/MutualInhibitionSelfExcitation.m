%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MutualInhibitionSelfExcitation.m
%A simple tri-stable loop for testing OLAC. 
%Four tunable parameters, 2 variables, three stables states. 
%To accompany "OptimalLeastActionControl_Top.m"

%Author: Daniel K. Wells (�) 2015
%Ver 1.0
%Email: dannykwells@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function out = MutualInhibitionSelfExcitation(t, x, P)
%This function is a tristable system representing a two genes which inhibit
%the other's expression, while promoting their own. 

%Note that everything in this function is vectorized. 

n1 = 4; k1 = 1; k2 = 1;  a1 = 1; a2 = 1; b1 = 1; b2 = 1;

%Baseline values for P: [.5 .5 .5 .5]; 


out1 = a1* x(1, :).^n1./(P(1)^n1 + x(1, :).^n1) + b1*(P(2)^n1)./(P(2)^n1 + x(2, :).^n1) - k1*x(1, :);
out2 = a2 *x(2, :).^n1./(P(3)^n1 + x(2, :).^n1) + b2*(P(4)^n1)./(P(4)^n1 + x(1, :).^n1) - k2*x(2, :);
out = [out1;out2];
end