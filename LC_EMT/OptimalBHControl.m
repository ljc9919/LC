function out=OptimalBHControl(Parameters) 


global  history

LC_Parameters = Parameters.LC_Parameters;

%Set
MaxIter = LC_Parameters.MaxIter; 
DiffMinChange = .01;%Set the DiffMinChange parameter of the optimization algorithm.
%Larger values should result in larger steps. 
DiffMaxChange = .5;%Set the DiffMaxChange parameters for the optimization algorithm. 
%Small values will result in smaller steps. 
ObjectiveLimit = LC_Parameters.ObjectiveLimit; 
TolCon = .01;%Set a tolerance on the constraint function. 
TolFunVal = 1e-4; %Set a tolerance on how little the function should change before the algorithm stops. 
TolXVal = 1e-4; %Set a tolerance on how little the X vector should change before the algorithm terminates. 
ConstrGradIncl = 'off';
ub = LC_Parameters.ub; 
lb = LC_Parameters.lb; 
InitParams =LC_Parameters.InitParams; 


%Set up the history.
history={}; history.ObjectiveValues={}; history.Params = {}; history.ConstraintViolations ={}; 
history.Parameters = Parameters; 


%Set up the top level optimization functions.
options = optimset('Algorithm', 'interior-point',...
'MaxIter',MaxIter,'Display', 'iter', 'MaxFunEvals', 1e10, ...
 'DiffMinChange', DiffMinChange,'DiffMaxChange', DiffMaxChange,  'ObjectiveLimit', ObjectiveLimit, 'TolCon', TolCon, ...
 'TolFun', TolFunVal, 'TolX', TolXVal, 'OutputFcn', @OutputFunction , 'GradConstr',ConstrGradIncl);

ConstraintFunction = Parameters.ConstraintFunction; 
ObjFunc = @(x)ObjectiveFunction(x, Parameters); 

[params, fval]=fmincon(ObjFunc, InitParams, [],[],[],[],lb, ub, ConstraintFunction, options); 
close all %Need to get rid of the window opened by history. 


out.ObjVal = fval; %What is the ultimate f value? 
out.Params = params; %What is the ultimate set of parameters?
out.History = history; %What is the history of the optimization? 



end




function fval = ObjectiveFunction(objvals, Parameters)
FP=Parameters.StableStates;
SP=Parameters.SaddleMatrix;
s_index=Parameters.IndexMatrix;

Function = Parameters.Function; 
% Jacobian = Parameters.Jacobian;
Functional = Parameters.Functional;


ID_Function = Parameters.ID_Function;

%Whole simulation constants
Dimension = size(FP, 1); 
NumFP = size(FP, 2);  


ContinuationParameters.in_params = Parameters.LC_Parameters.InitParams; 
ContinuationParameters.out_params = objvals; 
NewFP = zeros(Dimension, NumFP); 


%Calculate the identity of the fixed points. 
FP_IDs = zeros(1, NumFP); 

for i=1:NumFP
    ContinuationParameters.fpinit = FP(:, i); 
    ContinuationParameters.mode = i; 
    NewFPi=FixedPointContinuer(Function, ContinuationParameters); 
    NewFP(:, i) = NewFPi; 
    FP_IDs(i) = ID_Function(NewFPi); 
end

%Gauss Simulation
[alpha,mu,sigma] = GaussSimulate(objvals);
U_FP = zeros(1,NumFP);
for i=1:NumFP
    U_FP(i) = real(Ucal(alpha, mu, sigma, NewFP(:,i)));
end

%update saddle point
NewSP = cell(NumFP); 
U_SP = zeros(NumFP);
for i=1:NumFP
    for j=1:NumFP
        if i==j
            continue
        elseif isempty(SP{i,j})
            continue            
        else
%             NewSP{i,j} = hiosd(SP{i,j},objvals);
            NewSP{i,j} = hiosd(SP{i,j},alpha,mu,sigma,s_index(i,j));
            U_SP(i,j) = real(Ucal(alpha, mu, sigma, NewSP{i,j}));
            if isnan(U_SP(i,j))
%                 disp(objvals)
%                 disp(NewSP{i,j})
                U_SP(i,j) = max(U_FP(i)+1,U_FP(j)+1);
            elseif U_SP(i,j)==Inf
                U_SP(i,j) = max(U_FP(i)+1,U_FP(j)+1);
            end
        end
    end
end         
    

bhmat = zeros(NumFP);

for i=1:NumFP
    for j=1:NumFP
        if i==j
            BHVal = Inf;
%         elseif FP_IDs(i) == FP_IDs(j)
            %Then a bifurcation has occurred: these are the same spot. 
%             BHVal = 0;
        elseif isempty(SP{i,j})
            BHVal = Inf;
        else
            BHVal = U_SP(i,j) - U_FP(i);
        end
        bhmat(i,j) = BHVal;
    end
end
if ismember(inf,U_FP)
    bhmat = zeros(NumFP);
end
% disp(U_FP)
% disp(U_SP)
% disp(bhmat)
[fval, BHs, Occupancies] = Functional(bhmat, FP_IDs);
% [fval, alpha, Occupancies] = Functional(alpha, FP_IDs);


fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
ObjectiveValue = fval;
ObjectiveVariables = objvals;
FixedPoints =  NewFP;

% load occupancy.mat
% occupancy(end+1)=fval;
% save occupancy.mat occupancy
% 
% load time.mat
% time(end+1)=t;
% save time.mat time
% 
% save NewSP NewSP;

end

function stop  = OutputFunction(x, optimValues, state)
global history
stop = false; 


switch state 
    case 'init'
        hold on
    case 'iter'
        history.ObjectiveValues = [history.ObjectiveValues; optimValues.fval]; 
        history.Params =[history.Params,  x]; 
        history.ConstraintViolations = [history.ConstraintViolations, optimValues.constrviolation]; 
        
        
        %save(['history.mat'], 'history'); %If desired, the history can be
        %saved as you go.
        
        %if desired, the figure can be updated as you go. 
%         if history.Parameters.MakeFigures == 1
%         MakeFigure(history.Parameters, x); 
%         end
    case 'done'
        hold off
end

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





function [alpha,mu,sigma0] = GaussSimulate(P)
d=1;
cycle_index=500;

IM=load ('matrix.txt');
N = size(IM,1); %The dimension of the system.
Hilln=zeros(N,N); %Hill functiom
S=zeros(N,N); % threshold
AB=zeros(N,N); %the scale factors for the activation and inhibition
for i=1:N
    for j=1:N
        Hilln(i,j)=4;
        if IM(i,j)==1 && i~=j
            AB(i,j) = 2; %a
            if i>10 && j<11
            S(i,j) = .22; %sb
            else
            S(i,j) = 6; %sa
            end
        elseif IM(i,j)==1 && i==j
            AB(i,j) = 10; %aa
            if i>10 && j<11
            S(i,j) = .22; %sb
            else
            S(i,j) = 6; %sa
            end
        elseif IM(i,j)==-1
            AB(i,j) = 2; %b
            if i>10 && j<11
            S(i,j) = .22; %sb
            else
            S(i,j) = 6; %sa
            end
        end
    end
end


%%Solve odes from different initial values

p_num = 1;%%%%%%parallel
xx=zeros(N,cycle_index/p_num,p_num);

%%Solve odes from different initial values
k=ones(1,16);
% k(1)=P(1);k(4)=P(2);k(6)=P(3);k(8)=P(4);k(12)=P(5);k(13)=P(6);k(14)=P(7);k(15)=P(8);
k(1:9)=P(1:9);k(11:15)=P(10:14);
parfor p=1:p_num
    xxx = zeros(cycle_index/p_num,N);
    for i=1:(cycle_index/p_num)
        x0=unifrnd(0,20,[N,1]);
        [t,x]=ode15s(@(t,x)force2(t,x,k,S,IM,Hilln,AB),[0,20],x0);
        newx=x(end,:);
        x=inf*ones(1,N);
        nn=0;
        while norm( x(end,:)-newx(end,:) ,2 )>1e-3 && nn<=1e2
            x=newx;
            [t,newx]=ode15s(@(t,x)force2(t,x,k,S,IM,Hilln,AB),[0,1],x(end,:));
            nn=nn+1;
        end
        xxx(i,:)=newx(end,:);
    end
    xx(:,:,p) = xxx';
end
xx = reshape(xx,N,cycle_index);
xx = xx';
%%Finding the stable points
 for q=1:(cycle_index-1)
     for p=(q+1):cycle_index
         if norm(xx(q,:)-xx(p,:),'fro')<10^-1
             xx(p,:)=xx(q,:);
         end
     end
 end
stable_point=unique(xx(:,:),'rows');
n=zeros(1,2);
sigma=zeros(size(xx,1),size(xx,2)^2);
for i=1:size(stable_point,1)
    [m]=find(xx(:,2)==stable_point(i,2));
%     if length(m)>=1
%         disp(strcat(num2str(stable_point(i,:)),' repeat ',num2str(length(m)),' times',' the location in the row xx is' ,mat2str(m)))
%     end
    n(i,1)=m(1);
    n(i,2)=length(m);
    %%%calculate the covariance of each stable state
     sig=calculate_sigma(xx(m(1),:),N,d,k,S,IM,Hilln,AB);  
     for j=1:length(m)
         sigma(m(j),:)=sig;
     end
end

m=size(n,1);
if(m~=1)
    for i=1:m
        tran=n(i,:);
        flag=i;
        if( i~=m )
            for j=i+1:m
                if( xx(n(j,1),1)>xx(tran(1),1) )
                    tran=n(j,:);
                    flag=j;
                end
            end
        end
        n(flag,:)=n(i,:);
        n(i,:)=tran;
    end
end

index=size(n,1);  %% The number of the stable states
% alpha_w=zeros(index,1);
% w_pdf=zeros(index,1);
alpha=zeros(index,1);  %% The weight of the stable states
sigma0=cell(index,1);  %% The covariance of the Gaussian density function
mu=zeros(index,N);  %% The mean value of the Gaussian density function


for i=1:index
    %The mean value of each stable state
    mu(i,:)=xx(n(i,1),:); 
    %The covariance of each stable state
    sigma0{i}=reshape(sigma(n(i,1),:),N,N)';  
    %The weight of each stable state
    alpha(i)=n(i,2)/sum(n(:,2)); 
%     w_pdf(i) = mvncdf(zeros(1,N),15*ones(1,N),mu(i,:),sigma0{i});
%     alpha_w(i) = alpha(i)/w_pdf(i);
end
% alpha = alpha_w/sum(alpha_w);


end

function f=force2(t,x,k,S,IM,Hilln,AB)
%This is an implementation of solving the stable states of the ODEs
%%x:current location
%%a, aa, b, sa, sb, k, n: the parameters
%%mat:the matrix representing the relationship between genes
[num,m]=size(x);
H=zeros(num,m,num); 
for i=1:num
    for j=1:num 
        if IM(i,j)==1
            H(i,:,j)=Ha(x(i,:),S(i,j),AB(i,j),Hilln(i,j));
        elseif IM(i,j)==-1
            H(i,:,j)=Hr(x(i,:),S(i,j),AB(i,j),Hilln(i,j));
        end
    end
end
for i=1:num
    f(i,:)=sum(H(:,:,i))-k(i)*x(i,:);
end
end
function H=Hr(X,S,b,n)
H=b.*S.^n./(X.^n+S.^n);
end
function H=Ha(X,S,a,n)
H=a.*X.^n./(S.^n+X.^n);
end

function [sig]=calculate_sigma(xx,N,d,k,S,IM,Hilln,AB)

num = N;

jacobi_matrix = zeros(num,num);
for i = 1:num
    for j = 1:num
        if IM(i,j)==0
            J_ij=0;
        end
        if IM(i,j)==1
            J_ij =  AB(i,j)*deri_Ha(xx(i),S(i,j),Hilln(i,j));
        end
        if IM(i,j)==-1
            J_ij =  AB(i,j)*deri_Hi(xx(i),S(i,j),Hilln(i,j));
        end
        jacobi_matrix(i,j) = J_ij;
    end
end

Ajac = jacobi_matrix - k.*eye(num);

kk=N;
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
function deri_value = deri_Ha(x,s,n)
deri_value = s.^n.*n.*x.^(n-1)./(s.^n+x.^n).^2;
end

function deri_value = deri_Hi(x,s,n)
deri_value = -s.^n.*n.*x.^(n-1)./(s.^n+x.^n).^2;
end
   
function z=multivariate_normal_distribution(x,x0,sigma,n)
z=1/((2*pi)^(n/2)*det(sigma)^(1/2))*exp(1)^(-0.5*(x-x0)'*sigma^(-1)*(x-x0));
end

function out=Ucal(alpha, mu, sigma, x)
[index,num] = size(mu);
z = zeros(1,index);
for j=1:index
    z(j) = multivariate_normal_distribution(x, mu(j,:)',sigma{j},num);
end

H = z * alpha;
out = -log(H);

end



function out=hiosd(x,alpha,mu,sigma,index) %hiosd(x,objvals)

F=@(x)ngrad(x,alpha,mu,sigma);
H=@(x,v,l)hv(x,v,l,alpha,mu,sigma);

n = length(x);

optionsini = hiosdoptions(1);
optionsini.dt = 0.4;       
optionsini.outputX = 0; 
optionsini.outputp = 1;  
optionsini.outputd = 0;
optionsini.maxiter = 11; 
optionsini.k = 16;   
% rng(0);
optionsini.V = randn(n, optionsini.k);

[V] = hiosd_method1_ini(H, x, optionsini); 

options = hiosdoptions(1);
options.dt = 0.5;
options.betat = 1e1;       
options.betau = 0.001;
options.gamamax = 0.05;
options.gamamin = 1e-1;
options.outputp = 1;  
options.maxiter = 1000;  
options.epsf = 1e-50;
  
options.k = index;  
options.V = V(:, [1:index]);
% x0=x0-0.1*V(:,1);
out=hiosd_method1(F, H, x, options);

end

function out=hiosd_method1(F,H,x,options)

n = length(x);
innp = @(x,y) (x'*y);
nor = @(x) sqrt(x'*x);

options.V(isnan(options.V))=0;

if isfield(options,'k')
    if isfield(options,'V')
        if options.k ~= size(options.V,2)
            options.k = size(options.V,2);
%             disp('HiOSD: k is reset as the column - number of V.');
        end
    end
elseif isfield(options,'V')
    options.k = size(options.V,2);
else
    options.k = 1;
%     disp('HiOSD: k is reset as 1.');
end
k = options.k;

% V = zeros(n,k);
if isfield(options,'V')
    if rank(options.V) < k
        V = randn(n,k);
    else
        V = options.V;
    end
else
    V = randn(n,k);
end

for i = 1 : k
    V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i)/nor(V(:,i));
end

l=1e-6;

if isfield(options,'minl')
    minl=options.minl;
else
    minl=1e-6;
end

if isfield(options,'gamamax')
    gamamin=options.gamamax;
else
    gamamin=0.05;
end

% if isfield(options,'gamamin')
%     gamamin=options.gamamin;
% else
%     gamamin=1e-6;
% end

if isfield(options,'maxiter')
    maxiter=options.maxiter;
else
    maxiter=1e2;
end

dt=options.dt;
epsf=options.epsf;

alpha (1 : k, 1) = 0;

iter = 1;
f = F (x);
gp = zeros(n,1);
Dx = f - 2 * V * innp(V,f); 
Dx = Dx * dt;

DV = zeros(n,k);
DW = zeros(n,k);
U = zeros(n,k);
W = zeros(n,k);
Wp = zeros(n,k);
Vp = zeros(n,k);
gam (1 : k, 1) = 0;
DV(:,i) = V(:,i)*dt ;

for iter=1:maxiter
   
    tmpp = 2 * V * innp(V,f);
    g = f - tmpp;
    Dg = g - gp;
    gp = g;
    bta = abs ( innp(Dx,Dx) / innp(Dx,Dg) );
    bta = min ( bta, options.betat);
    bta = max ( bta, options.betau);
    
    xp = x;
    x = xp + bta * g;
    l = max(l / (1 + bta),minl);
    Dx = x - xp;
    
    for i = 1 : k
        U(:,i) = H (x, V(:,i), l);
        alpha(i) = innp (V(:,i), U(:,i));
        W(:,i) =  ( U(:,i) - alpha(i) * V(:,i) ) ;
        W(:,i) = W(:,i) - 2 * V(:,1:i-1) * innp( V(:,1:i-1), W(:,i) ) ;
        DW(:,i) = W(:,i) - Wp(:,i);
        Wp(:,i) = W(:,i);
        gam(i) = abs ( innp(DV(:,i),DV(:,i)) / innp(DV(:,i),DW(:,i)) );
%         gam(i) = max ( gam(i), gamamax);
        gam(i) = max ( gam(i), gamamin);
    end
    
    for i = 1 : k
        Vp(:,i) = V(:,i);
        V(:,i) = V(:,i) - gam(i) * W(:,i);
        V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
        V(:,i) = V(:,i)/nor(V(:,i));
        DV(:,i) = V(:,i) - Vp(:,i);
    end
    
    f = F(x);
    res = norm(f);
    
%     if res < epsf || res > 1e+20
    if res < epsf
        break;
    end
%     if mod(iter,options.outputp)==0
%     disp([num2str(iter) '   ' num2str(res)]);
%     disp(x);

% output = struct('x', x, 'V', V, 'it', iter,'res',res);
% output = struct('x', x, 'it', iter,'res',res);
end
out = x;
% disp(x)
% disp(res)
% disp(V)
end

% function F = ngrad(x,P)
% n1 = 4; k1 = 1; k2 = 1;  a1 = 1; a2 = 1; b1 = 1; b2 = 1;
% F = [a1* x(1, :).^n1./(P(1)^n1 + x(1, :).^n1) + b1*(P(2)^n1)./(P(2)^n1 + x(2, :).^n1) - k1*x(1, :);...
%      a2 *x(2, :).^n1./(P(3)^n1 + x(2, :).^n1) + b2*(P(4)^n1)./(P(4)^n1 + x(1, :).^n1) - k2*x(2, :)];
% end
% 
% function hv = hv(x,v,l,P)
% F = @(x)ngrad(x,P);
% hv= -(F(x+l*v) - F(x-l*v)) ./ (2*l);
% end

function F=ngrad(x,alpha,mu,sigma)
[n,dim] = size(mu);
grad = zeros(dim,n);

for i=1:n
    grad(:,i) = gaussgrad(x,mu(i,:)',sigma{i},dim);
end

F = grad * alpha; 
end

function g=gaussgrad(x,mu,sigma,n)
dim = length(mu);
g = zeros(dim,1);
f = 1/((2*pi)^(n/2)*det(sigma)^(1/2))*exp(1)^(-0.5*(x-mu)'*sigma^(-1)*(x-mu));
sigma_inv = inv(sigma);
for i=1:dim
    gp = 0;
    for j=1:dim
        gp = gp + (x(j)-mu(j))*sigma_inv(i,j);
    end
    g(i) = -f * gp * 0.5;
end
% g = [f*(-1/det(sigma)*(sigma(2,2)*(x(1)-mu(1))-sigma(1,2)*(x(2)-mu(2))));
%      f*(-1/det(sigma)*(sigma(1,1)*(x(2)-mu(2))-sigma(1,2)*(x(1)-mu(1))))];
end

function hv = hv(x,v,l,alpha,mu,sigma)
F = @(x)ngrad(x,alpha,mu,sigma);
hv= -(F(x+l*v) - F(x-l*v)) ./ (2*l);
end


function [V, alpha, iter] = hiosd_method1_ini(H,x,options)

n=length(x);
innp = @(x,y)(x'*y);
nor = @(x)sqrt(x'*x);

if isfield(options,'k')
    if isfield(options,'V')
        if options.k ~= size(options.V,2)
            options.k = size(options.V,2);
%             disp('HiOSD: k is reset as the column - number of V.');
        end
    end
elseif isfield(options,'V')
    options.k = size(options.V,2);
else
    options.k = 1;
%     disp('HiOSD: k is reset as 1.');
end
k = options.k;

V = zeros(n,k);

if isfield(options,'V')
    if rank(options.V) < k
        V = randn(n,k);
    else
        V = options.V;
    end
else
    V = randn(n,k);
end

for i = 1 : k
    V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
    V(:,i) = V(:,i)/nor(V(:,i));
end

l=1e-6;

if isfield(options, 'minl')
    minl=options.minl;
else
    minl=1e-6;
end


if isfield(options,'gamamin')
    gamamin=options.gamamin;
else
    gamamin=1e-6;
end

if isfield(options,'maxiter')
    maxiter=options.maxiter;
else
    maxiter=1e3;
end

dt=options.dt;

iter = 1;
DV = zeros(n,k);
DW = zeros(n,k);
U = zeros(n,k);
W = zeros(n,k);
Wp = zeros(n,k);
Vp = zeros(n,k);
gam (1 : k, 1) = 0;
DV(:,i) = V(:,i) * dt;

for iter=1:maxiter
    for i = 1 : k
        U(:,i) = H (x, V(:,i), l);
        alpha(i) = innp (V(:,i), U(:,i));
        W(:,i) =  ( U(:,i) - alpha(i) * V(:,i) ) ;
        W(:,i) = W(:,i) - 2 * V(:,1:i-1) * innp( V(:,1:i-1), W(:,i) ) ;
        DW(:,i) = W(:,i) - Wp(:,i);
        Wp(:,i) = W(:,i);
        gam(i) = abs ( innp(DV(:,i),DV(:,i)) / innp(DV(:,i),DW(:,i)) );
        gam(i) = max ( gam(i), gamamin );
    end
    
     for i = 1 : k
        Vp(:,i) = V(:,i);
        V(:,i) = V(:,i) - gam(i) * W(:,i);
        V(:,i) = V(:,i) - V(:,1:i-1) * innp(V(:,1:i-1), V(:,i));
        V(:,i) = V(:,i)/nor(V(:,i));
        DV(:,i) = V(:,i) - Vp(:,i);
     end
    
%     if mod(iter, 1)==0
%         disp(V);
%     end
    iter = iter+1;
%     disp(alpha')
end

end

function options = hiosdoptions(x)
%% HiOSD options setting
% x       initial value of x, column vector n*1
% V       initial value of V; n*k
% innpfunc can compute multi - one column inner product as a column vector
% k       index of the sddle point
% dt      initial step size 
%         if dt = 0, then the iteration of x is stopped
% l       initial dimer length
% minl    minimal dimer length
% epsf    iteration stops if |F(x)| is smaller than this
% maxiter max iteration step number
% betat   max beta
% betau   min beta
% gamamax max gamma
% gamamin min gamma
% outputX store iteration message every outputX steps
% outputp print iteration message every outputp steps
% outputd draw iteration message every outputp steps
% draws   draw sentence
% tau     max |x-xp| to control step size
% precd   preconditioner for LOBPCG method
% adding  additional sentence every iteration

options=struct('k',1,'dt',0.1,'l',1e-6,'minl',1e-6,...
    'epsf',1e-10,'maxiter',1e3,...
    'betat',max(x,1),'betau',min(1,1/x),'gammamax',max(x,1),'gammamin',min(1,1/x),...
    'outputX',1, 'outputp',1, 'outputd', 1, 'draws', '',...
    'tau',0.5,'precd',@(p)p);

if x==Inf
    options.betat=Inf;
    options.betau=0.01;
    options.gammammax=Inf;
    options.gammammin=1;
end
end


