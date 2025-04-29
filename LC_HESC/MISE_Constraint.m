function [C_ineq, C_eq, Cineq_grad, Ceq_grad] = MISE_Constraint(x, ConstraintParams)
%x should be the tunable parameters of the function. 
%This function includes all constraints on the point x, including all
%bifurcation conditions, etc. 

%Here we use this function to ensure that specified stable fixed points
%stay in existence. 

Function = ConstraintParams.Function; %This should be the main function, and it should 
%include a time component as well. 

Jacobian = ConstraintParams.Jacobian; 
ID_Function = ConstraintParams.ID_Function;

original_fp = ConstraintParams.OriginalFP; 
original_params = ConstraintParams.OriginalParams; 
% L = ConstraintParams.L;


BifurcationVector = [1 1 1];  %A 1,0 array of which fp should be preserved.
%1: Preserve the specified state
%2. The specified state can undergo a bifurcation. 
%This vector can be modified if it is acceptable for certain states to be
%lost to bifurcation. I


l_Bifurcation = sum(BifurcationVector); 
C_ineq = zeros(1, l_Bifurcation+1); 
tol =ConstraintParams.Tolerance ; %The largest real eigenvalue should thus not be larger than -.05; 



 
inparams.in_params = original_params;
inparams.out_params = x; 
 


for i=1:l_Bifurcation
    inparams.fpinit = original_fp(:, i); 
    C_ineq(i) = Adaptive_FixedPointContinuer_Bifurcations(Function,Jacobian, inparams) + tol;     
end


ContinuationParameters.in_params = original_params; 
ContinuationParameters.out_params = x; 

% lambda=sum(abs(x-original_params));
% C_ineq(l_Bifurcation+1)=lambda-L;

new_fp = zeros(size(original_fp));
%Calculate the identity of the fixed points. 
NumFP = size(new_fp,2);
fp_IDs = zeros(1,NumFP); 

for i=1:NumFP
    ContinuationParameters.fpinit = original_fp(:, i); 
    ContinuationParameters.mode = i; 
    NewFPi=FixedPointContinuer(Function, ContinuationParameters); 
    new_fp(:, i) = NewFPi; 
    fp_IDs(i) = ID_Function(NewFPi); 
end
n_stable = length(unique(fp_IDs));
C_eq = n_stable - NumFP; 

% C_ineq_p = zeros(1,NumFP);
% for i=1:NumFP
%     sig = fliplr(calculate_sigma(new_fp(:,i), x, length(new_fp(:,i)), 0.05));
%     sigma = reshape(sig,2,2);
%     C_ineq_p(i) = -prod(eig(sigma));
% end
% 
% C_ineq = cat(2, C_ineq, C_ineq_p);

Cineq_grad = []; 
Ceq_grad = []; 


end


function [maxRealEigenPart] = Adaptive_FixedPointContinuer_Bifurcations(func,J, params)
%This is a straight forward method to calculate the largest real part the
%the eigenvalues associated with a particular stable state as that state is
%continued. Since this value changes very suddenly as a bifurcation is
%appraoched, an adaptive scheme had to be implemented. 

format long; 

%What is the fixed point that we want to continue?
fpinit = params.fpinit; 
originalfp = fpinit; 
 
%What is the initial set of parameters that corresponds to this fixed point?
in_params = params.in_params; 
num_params = size(in_params, 2); 

%What set of parameters are we going to? 
out_params=params.out_params; 

%How many steps should the continuation take? 
N1 = 31; %Number of steps to do the initial continuation
N2 = 31; %Number of steps to do teh adaptive continuation of. 
 

%Calculate the line for the homotopy method. 
param_line = zeros(num_params, N1); 
for i=1:num_params
    param_line(i, :) = linspace(in_params(i), out_params(i), N1); 
end



%This N should be relativly small. 
norm_mat = zeros(16, N1+1); 
norm_mat(:, 1) = fpinit;
for i=1:N1
    params_vals=param_line(:, i); 
    f = @(t, x)func(t, x, params_vals);  
    %Calculate the new fixed point
    [T, Y]=ode15s(f, [0, 1e3], fpinit);
    fpinit_test=Y(end, :)';
    norm_mat(:,i+1) = fpinit_test; %Save the change in norm between fixed point. 
    %Those places where the stable state is changing the most are where a
    %bifurcation could be occurring. 
    fpinit = fpinit_test; 
end
% 

%This is an adaptive scheme to calculate the largest real part of the
%jacobian of the system where the stable states are changing the most. 
%This is essentially a numerical variational problem. 
%dparams = out_params - in_params; %These are the changes in each fixed point. 
dnorm = sum(abs((norm_mat(:, 2:end) - norm_mat(:, 1:end-1))));  


w = sqrt(1 + 1e10*dnorm.^2); 
alphak = [0 cumsum(w)/sum(w)];
remesh= interp1(alphak, 0:N1, linspace(0, 1, N2));
remesh = remesh/max(remesh); 


newparams = zeros(size(in_params, 2), N2); 
for i=1:size(in_params, 2)
    newparams(i, :) = spline(linspace(0, 1, N1),param_line(i, :), remesh); 
end


Eigen = -1; 
%EigMat = zeros(1, N2);
fpinit = originalfp; 
for i=1:N2
    params_vals=newparams(:, i); 
    

    f = @(t, x)func(t, x, params_vals);  
    jtemp = @(x)J(1, x, params_vals); 

    %Calculate the new fixed point
    [T, Y]=ode15s(f, [0, 1e3], fpinit);
    fpinit_test=Y(end, :)';
    warning off
    Eigentest = max(real(eigs(jtemp(fpinit)))); 
    Eigen = max(Eigen, Eigentest); 
    %EigMat(i) = Eigen; 
    fpinit = fpinit_test; 
    
end
% 



maxRealEigenPart = Eigen;  
end

function [out] = FixedPointContinuer(func, params)
%This is a straight forward method to do a possibly high dimensional
%homotopy method of a stable fixed point. It checks the determinant to make
%sure that a bifurcation is not occurring, but aside from that doesn't do
%anything fancy to work around it. 

format long; 


%What is the fixed point that we want to continue?
fpinit = params.fpinit; 

%What is the initial set of parameters that corresponds to this fixed point?
in_params = params.in_params; 
num_params = size(in_params, 2); 

%What set of parameters are we going to? 
out_params=params.out_params; 

%How many steps should the continuation take? 
N=50;  


%Calculate the line for the homotopy method. 
param_line = zeros(num_params, N); 
for i=1:num_params
    param_line(i, :) = linspace(in_params(i), out_params(i), N); 
end

  
for i=1:N
    params_vals=param_line(:, i); 

    f = @(t, x)func(t, x, params_vals); 

    %Calculate the new fixed point
    [T, Y]=ode15s(f, [0, 1e3], fpinit);
    fpinit_test=Y(end, :)';

    fpinit = fpinit_test; 
end

out=fpinit; 
end

function [sig]=calculate_sigma(xx,par,kk,d)
% alphaP=par(1);alphaA=par(2);betaP=par(3);betaA=par(4);rP=par(5);rA=par(6);h=par(7);uP=par(8);
% uA=par(9);K=par(10);betaother=par(11);r0=par(12);epsilon=par(13);
% global r1 r2 net

syms x1 x2
% for i=1:size(net,1)
%     P(i)=sym(['P' num2str(i)]);
% end
% 
% syms A
% for i=1:size(net,2)
%     A(i)=sym(['A' num2str(i)]);
% end

for i=1:kk
    mm(i)=xx(i);
end

% for i=1:size(net,2)
%     nn(i)=xx(i+size(net,1));
% end
   
% Ajac=jacobian([alphaP*P'-betaP*P'.*P'-betaother*P'*sum(P)+P'.*r1*A'./(1+h*r1*A')+uP;...
%      alphaA*A'-betaA*A'.*A'-betaother*A'*sum(A)-K*A'+A'.*r2'*P'./(1+h*r2'*P')+uA],...
%      [P,A]);
Ajac=jacobian([(par(1) * x1^4)/(par(4)^4 + x1^4) + (par(2) * par(4)^4)/(par(4)^4 + x2^4) - x1;...
               (par(1) * x2^4)/(par(4)^4 + x2^4) + (par(2) * par(4)^4)/(par(4)^4 + x1^4) - x2;],[x1,x2]);
% Ajac=jacobian([12000 * hill(x1,2,6000,3) * hill(x2,0.1,9000,3) * hill(x3,0.5,2000,3) - 1.04 * x1;
%                9000 * hill(7600,2,25000,3) * hill(x2,1.25,13000,3) * hill(x3,0.5,4000,3) * hill(x1,0.2,6000,5) - 0.5 * x2;
%                400 * hill(x3,2,6500,4) * hill(x2,1.4,18000,3) - 0.1 * x3;],[x1,x2,x3]);

Ajac=double(subs(Ajac,[x1,x2],[mm]));

% A*sigma+sigma*A'+2D

P=zeros(kk^2,kk^2);  %coefficient matrix

%%the initial of coeffiicient matrix
for i=0:(kk-1)
    P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)=P(i*kk+1:i*kk+kk,i*kk+1:i*kk+kk)+Ajac;
end

for m=0:kk-1
    for i=1:kk
        for j=1:kk
            P(m*kk+i,(j-1)*kk+i)=P(m*kk+i,(j-1)*kk+i)+Ajac(m+1,j);
        end
    end
end

B=zeros(kk^2,1);
for i=1:kk
    B((i-1)*kk+i)=-2*d;
end


sig=P\B;

end