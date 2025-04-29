%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MISE_Jacobian.m
%The Jacobian of the dynamical system included. 
%Note that this function is vectorized. In general, any Jacobian given to
%AMAM.m must be vectorized in the third dimension to be usable. 


%Author: Daniel K. Wells (©) 2015
%Ver 1.0
%Email: dannykwells@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function out = MISE_Jacobian(t, xin, P)
if ndims(xin)<3
    in = zeros(1, size(xin, 1), size(xin, 2)); 
    in(1, :, :) = xin; 
    xin=in; 
end
tsteps = size(xin, 3); %If there is a constant or a zero in the jacobian...
format long

x1 = xin(:, 1, :); 
x2 = xin(:, 2, :); 

n = 4; k1 = 1; k2 = 1;  a1 = 1; a2 = 1; b1 = 1; b2 = 1;

out = [x1.^(-1).*(P(1).^n+x1.^n).^(-2).*(a1.*n.*P(1).^n.*x1.^n+(-1).*k1.*x1.*(P(1).^n+x1.^n).^2),...
    (-1).*b1.*n.*P(2).^n.*x2.^((-1)+n).*(P(2).^n+x2.^n).^(-2);...
    (-1).*b2.*n.*P(4).^n.*x1.^((-1)+n).*(P(4).^n+x1.^n).^(-2),...
    x2.^(-1).*(P(3).^n+x2.^n).^(-2).*(a2.*n.*P(3).^n.*x2.^n+(-1).*k2.*x2.*(P(3).^n+x2.^n).^2)];


end