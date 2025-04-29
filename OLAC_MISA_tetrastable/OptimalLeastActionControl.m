%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is our implementation of optimal least action control. 
%The license of the top level function (OptimalLeastActionControl_Top.m) holds for this function and all
%subfunctions. 

%Author: Daniel K. Wells (�) 2015
%Ver 1.0
%Email: dannykwells@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=OptimalLeastActionControl(Parameters) 


global  history

OLAC_Parameters = Parameters.OLAC_Parameters;

%Set
MaxIter = OLAC_Parameters.MaxIter; 
DiffMinChange = .05;%Set the DiffMinChange parameter of the optimization algorithm.
%Larger values should result in larger steps. 
DiffMaxChange = .2;%Set the DiffMaxChange parameters for the optimization algorithm. 
%Small values will result in smaller steps. 
ObjectiveLimit = OLAC_Parameters.ObjectiveLimit; 
TolCon = .01;%Set a tolerance on the constraint function. 
TolFunVal = 1e-3; %Set a tolerance on how little the function should change before the algorithm stops. 
TolXVal = 1e-5; %Set a tolerance on how little the X vector should change before the algorithm terminates. 
ConstrGradIncl = 'off';
ub = OLAC_Parameters.ub; 
lb = OLAC_Parameters.lb; 
InitParams =OLAC_Parameters.InitParams; 


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

Function = Parameters.Function; 
Jacobian = Parameters.Jacobian;
Functional = Parameters.Functional;


ID_Function = Parameters.ID_Function;

%Whole simulation constants
Dimension = size(FP, 1); 
NumFP = size(FP, 2);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AMAM Parameters
AMAM_Params=AMAM_params(Parameters);

%OLAC relies on the Adaptive Minimum Action Method to find the minimum of
%the action functional. These parameters should NOT need extensive tuning
%for different applications. 


ContinuationParameters.in_params = Parameters.OLAC_Parameters.InitParams; 
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




Func = @(x)Function(1, x, objvals); 
dFunc = @(x)Jacobian(1, x, objvals); 


actionmat = zeros(NumFP);

for i=1:NumFP 
    for j=1:NumFP
        if i==j
            ActionVal = Inf; 
        elseif FP_IDs(i) == FP_IDs(j)
            %Then a bifurcation has occurred: these are the same spot. 
            ActionVal = 0; 
        else
            Spot1 = i; Spot2 = j;

            AMAM_Params.init = NewFP(:, [Spot1,Spot2]);

            [ActionVal, Path] = AMAM(AMAM_Params, Func, dFunc);  

        end
        actionmat(i,j) =ActionVal;
    end
end

[fval, Actions, Occupancies] = Functional(actionmat, FP_IDs); 

fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
ObjectiveValue = fval;
ObjectiveVariables = objvals;
FixedPoints =  NewFP;




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
        if history.Parameters.MakeFigures == 1
        MakeFigure(history.Parameters, x); 
        end
    case 'done'
        hold off
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FixedPointContinuer.m
%Given a dynamical system and an start/end parameter set, continue a specified
%stable state from the starting point to the ending one. 
%We rely on a simple homotopy method to accomplish this. More advanced
%methods are of course possible, which is why it is recommended to include
%a bifurcation monitoring function. 
%Without, the function will check the determinant of the
%jacobian of the system for its distance to 0. 

%Author: Daniel K. Wells (�)
%Date: 11/11/2014
%Ver 1.0
%Email: dannykwells@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AMAM.m
%This is an implementation of the adaptive minimum action method developed
%by Zhou et al. in 
%Zhou X, Ren W, E W (2008) An adaptive minimum action method for the study
%of rare events. J. Chem Phys 128:104111. 
%It is to accompany the software suite with top level
%"OptimalLeastActionControl_Top.". All license information in that file
%pertains to this file as well. 
%This is sparsely documented. A more thorough description of the method can
%be found in the above paper. 
%For use with Optimal Least Action Control, nothing in this code should be
%altered or changed, unless a bug is found. 

%Author: Daniel K. Wells (�)
%Date: 11/11/2014
%Ver 1.0
%Email: dannykwells@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [fval, phi] = AMAM(params, funcs, dfuncs)
format long
 
%The Jacobian should be such that it can take in a 2-d array of values
%through time, and return the 3-d time Jacobian. Time should be on the
%third axis. 


xffp = params.init;
%The fixed points should be given in columns. 
tmax = params.TMax;
n = params.N; 
%Note that n is the number of points to use, not the dimension.
m = params.Dimension;
c = params.c;
maxiter = params.MaxIter;
kmax = params.K;
remesh = params.Remesh;
qmin = params.q;
phi = params.PhiInit;

phil = xffp(:,1); 
phir = xffp(:,2);

T = -tmax/2:tmax/n:tmax/2;
delt = tmax/n*ones(m, n);       
k = 0;  

%If an initial path is not set, use a straight line. 
if isempty(params.PhiInit) 
    phi = pathfinder(phil, phir, T);
end

while k < kmax 

    %ub is alwsy set to be inf. 
     ub = inf*ones(m,n-1);
     lb = -inf*ones(m, n-1); 
      

     if params.ObjGrad==1     
         options = optimset('Algorithm', 'interior-point','Hessian',{'lbfgs',5},...
        'MaxIter',maxiter,  'MaxFunEvals', 3000000,...
    'AlwaysHonorConstraints', 'none','Display', 'off', 'GradObj','on');
        func = @(x)Sg(x, delt,funcs,dfuncs,phir, phil);
        
     else
         options = optimset('Algorithm', 'interior-point','Hessian',{'lbfgs',5},...
        'MaxIter',maxiter,  'MaxFunEvals', 3000000,...
       ...% 'TolX', 1e-10,'TolCon', 1e-10, 'TolFun', 1e-10,
    'AlwaysHonorConstraints', 'none','Display', 'off');
        func = @(x)S(x, delt,funcs,phir, phil); 
       
     end

    phi = phi(:, 2:end-1);



    [x2, fval] = fmincon(func, phi, [],[],[],[],...
                         lb,ub,[],options);

     x2mod = [phil, x2, phir];
     w = monitor(x2mod, delt, c);
     q = qfunc(w, delt);
     if remesh == 1 && q > qmin
         [T, phi, delt] = gridremesh(n,m,w,delt,T, x2mod);
         else
         phi = x2mod;
     end
     k = k+1  ; 
end
phi = phi(:, 2:end-1); 
         
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [out, g] = Sg(phi, delt, func,dfunc, phirp, philp)
%Approximate the action functional according to the midpoint
%rule. The time derivative of phi is taken inside the function.

%Create the relevant phi vectors
phi = [philp,phi, phirp];
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phihalf = (phir + phil)/2;




%Perform the trapezoidal approximation. 
summand = sum(((phir - phil)./(delt) - func(phihalf)).^2);
%summand = sum(summand.^2);             
out = .5*sum(summand.*delt(1, :));

g = gradS(phi(:,2:end-1) , delt, func, dfunc, phirp, philp); 
end

function  out = S(phi, delt, func, phirp, philp)
%Approximate the action functional according to the midpoint
%rule. The time derivative of phi is taken inside the function.

%Create the relevant phi vectors
phi = [philp,phi, phirp];
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phihalf = (phir + phil)/2;




%Perform the trapezoidal approximation. 
summand = sum(((phir - phil)./(delt) - func(phihalf)).^2);
%summand = sum(summand.^2);             
out = .5*sum(summand.*delt(1, :));



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = monitor(phi, delt, c)
%The monitor function
format long
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phit = (phir - phil)./delt;
derivmagnitude = sum(phit.^2);
w = sqrt(1 + c*derivmagnitude);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tnew, phi, delt] = gridremesh(n,m,w,delt,T,x2mod)
%The adaptive remeshing procedure, as in E et al. 
format short
%Calculated the integral of w using the midpoint rule. 
intw = sum(w.*delt(1, :));
%Calculate alpha in accordance with the paper.
alphak = [0,cumsum((w.*delt(1, :))/intw)];
%The uniform alpha. 
alphaold = 0:1/n:1;
%This is a linear interpolation of T as a function of alpha from the
%stretched alphak to the uniform alpha. 
Tnew = interp1(alphak, T, alphaold, 'linear');
%interp1 was giving a NaN in this position for reasons unknown, since it
%should just give the end value back. This is an ragged solution to the
%problem. 
Tnew(end) = -T(1);
%The matrix of delt. 
delt = ones(m, 1)*(Tnew(2:end) - Tnew(1:end-1));
%Find the new phi by a cubic interpolation over the remeshed T.
phi = spline(T, x2mod, Tnew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = qfunc(w, delt)
%The calculation of q, in accordance with E et al. 
       maxx = max(w.*delt(1, :)); minx = min(w.*delt(1, :));
       q = maxx/minx;
       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path = pathfinder(start, ending, time)
%given to starting points, creates a matrix of a linear transit between
%them. 

dims = size(start);
dim = dims(1);
time = time + time(end);

% [xstart, ystart] = meshgrid(start);
% [xend, yend] = mehsgrid(end);
% [xt, yt]= meshgrid(time);
% 
% path

xstart = start*ones(1, length(time));
diff = (ending- start)*ones(1, length(time));
timep = (time/(max(time)))';
times = (timep*ones(1, dim))';
path = xstart + diff.*times;


end

function out = gradS(phi, delt, func, dfunc, phirp, philp)

%Add down, not across, for the gradient. 
phi = [philp,phi, phirp];
phir = phi(:,2:end);
phil = phi(:,1:end-1);
phihalf = (phir + phil)/2;
%These go across all dimensions and time.
sz = size(phihalf); 
% out = zeros(sz(1), sz(2)-1); 
% ndim = sz(1); 
%For now, we will have to use a single loop
T0 = ((phir - phil)./(delt) - func(phihalf)); 

T0p = reshape(T0, [1,sz(1), sz(2)]); 
deltp = reshape(delt,  [1, sz(1),  sz(2)]); 
phijacob = dfunc(phihalf);

out = zeros(1, sz(1), sz(2)-1); 
for i=1:sz(1); 
    a1 = T0p(1, i, 1:end-1); 
    a2 = T0p(1, i, 2:end); 
    
    step1 = a1.*(1./(deltp(1, i, 1:end-1)) - 1/2*phijacob(i,i,1:end-1)).*deltp(1, i, 1:end-1);
    step2 = a2.*(-1./(deltp(1, i, 2:end)) - 1/2*phijacob(i,i,2:end)).*deltp(1, i, 2:end);
    vals = 1:sz(1); 
    vals(i) = []; 
    step3 = -1/2*(sum( T0p(1, vals, 1:end-1).*permute(phijacob(vals,i, 1:end-1), [2, 1, 3]).*deltp(1, vals, 1:end-1), 2) + sum( T0p(1, vals, 2:end).*permute(phijacob(vals,i, 2:end), [2, 1, 3]).*deltp(1, vals, 2:end), 2)); 

    out(1, i, :) = step1 + step2 + step3; 
    
end

out = squeeze(out); 


end


function AMAM_Params=AMAM_params(Parameters)
%These parameter can be modified if desired, but they are not included for
%modification above. 
AMAM_Params={}; 
AMAM_Params.N= 15; %The number of points in the minimum action path. 
%Larger is more accurate, but slower. 
AMAM_Params.TMax = 55; %The time range of the minimum action path. 
%Larger values can be more accurate, but can also lead to instabilities. 
AMAM_Params.Dimension = size(Parameters.StableStates, 1); %The dimension of the system. 
AMAM_Params.c = 1e10; %The remeshing parameter, c. Larger values give greater remeshing. 
AMAM_Params.MaxIter = 300; %The number of optimization steps to run between remeshing steps. 
AMAM_Params.K= 2; %The number of total remeshing steps to do; 
AMAM_Params.Remesh = 1; %Turn the remesh on or off. If off, the algorithm will be default perform K*MaxIter itereations. 
AMAM_Params.q = 3; %The q parameter from the original paper. Is a measure of path smoothness. 
AMAM_Params.ub = []; %If desired, constraints can be set of the path as well. This is not done in our manuscript. 
AMAM_Params.lb = []; %As above, lower bounds can be set on the path as well. This is also not done at all in our manuscript. 
AMAM_Params.PhiInit =[]; %The default intial path between two states is a straight line. This can be modified here. 
AMAM_Params.ObjGrad = Parameters.JacobianPresent;%We require a Jacobian in this implementation. 
    
end

function MakeFigure(Parameters, out)



%%%%%%
%Plot the final picture
FP = Parameters.StableStates; 
Params = Parameters.OLAC_Parameters.InitParams; 
Function = Parameters.Function; 
Jacobian = Parameters.Jacobian;
Func = @(x)Function(1, x, out); 
dFunc = @(x)Jacobian(1, x, out); 
szFP= size(FP, 2); 
AMAM_Params=AMAM_params(Parameters);
colors = {[203/255, 50/255, 219/255] , ...
    [35/255,222/255, 52/255], [26/255, 157/255, 201/255], [227/255 124/255 47/255]};
place = 1; 
StorageCell = cell(2, szFP^2-szFP); 
ActionMat  =zeros(szFP); 
FP_Cell = cell(1, szFP);
close all
for ii=1:size(FP, 2); 
    for jj=1:size(FP, 2); 
        if ii~=jj
            
        ContinuerParams.fpinit = FP(:, ii); 
        ContinuerParams.in_params = Params; 
        ContinuerParams.out_params = out;
        FPii = FixedPointContinuer(Function, ContinuerParams);
        
        
        ContinuerParams.fpinit = FP(:, jj); 
        FPjj = FixedPointContinuer(Function, ContinuerParams);
            
        FP_Cell{ii} = FPii; 
            
        AMAM_Params.init = [FPii, FPjj]; 
        AMAM_Params.TMax = 55;
        AMAM_Params.N = 15; 
        AMAM_Params.Remesh = 1; 
        AMAM_Params.K= 3;
        [ActionVal, Path] = AMAM(AMAM_Params, Func, dFunc);
        StorageCell{1,place} = Path; 
        StorageCell{2, place} = ActionVal; 
        StorageCell{3, place}=[ii,jj];
        place = place +1;
        
        ActionMat(ii,jj) = ActionVal; 
        end
        

    end   
end

Num = place-1;
NewStorageCell = StorageCell; 
Places2Delete = [];
for ii=1:Num
    for jj=1:Num          
        if (ActionMat(StorageCell{3, ii}(1),StorageCell{3,ii}(2) ) + ActionMat(StorageCell{3, jj}(1),StorageCell{3,jj}(2) )) < ActionMat(StorageCell{3, ii}(1),StorageCell{3,jj}(2) )
            for kk=1:Num
                if StorageCell{3, kk} == [StorageCell{3, ii}(1) StorageCell{3, jj}(2)]
                    Places2Delete = [kk Places2Delete];  
                end
            end
        end
    end
end
for ii=1:Num
    if StorageCell{2, ii} == 0
        Places2Delete = [ii Places2Delete]; 
    end
end

NewStorageCell(:, unique(Places2Delete))=[];
Num = size(NewStorageCell, 2); 

X = -.1:.2:2.1; Y = -.1:.2:2.1;
[i,j] = meshgrid(X,Y);
cf = @(x, y)Func([x;y]);
arrows = cf(i(:)', j(:)');

figure(1)
h(1) = quiver(i,j,reshape(arrows(1, :), 12, 12), reshape(arrows(2, :),12,12), 'LineWidth', 1.2,...
    'Color',[72/255, 92/255, 99/255] ); hold on;
axis( [-0.1 2.2 -0.1 2.2])

for ii=1:Num
    plot(NewStorageCell{1,ii}(1, :), NewStorageCell{1,ii}(2, :), 'LineWidth', 2.5, 'Color', colors{ii});
end

for i=1:length(FP_Cell)
    
rectangle('Position',[FP_Cell{i}(1)-.05,FP_Cell{i}(2)-.05,.1,.1],...
    'Curvature',[1,1],...
    'FaceColor','black'), hold on
daspect([1,1,1])
end
    
    
set(gca, 'FontSize', 14)
xlabel('x_1', 'FontSize', 14); 
ylabel('x_2', 'FontSize', 14);
pause(1); 





end




