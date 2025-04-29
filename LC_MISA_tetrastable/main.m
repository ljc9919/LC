%Define the parameter struct on which the algorithm will operate. 
Parameters = {}; 



%1. Define the (controllable dynamical system). 
%This should be a function that takes in three parameters, in this order: 
%(t, x, P), where t is the time (for non-autonomous dynamical systems), x
%is the state vector, and P is the parameter vector for the system. 


Function  = @(t, x, P) MutualInhibitionSelfExcitation(t, x, P);
Parameters.Function = Function;



%2. The states of the system must be previously identified to run LC. These
%should be given as an n-row by m-column array, where n is the number of
%dimensions of the dynamical system under consideration, and m is the
%number of stable states. 
%should be inputted here:

StableStates = load('MISE_StableStates4');
StableStateNames = fieldnames(StableStates); 
StableStates = StableStates.(StableStateNames{1}); 
Parameters.StableStates = StableStates; 


SaddleMatrix = load('MISE_SaddleMatrix_Gauss4');
SaddleMatrixNames = fieldnames(SaddleMatrix); 
SaddleMatrix = SaddleMatrix.(SaddleMatrixNames{1}); 
Parameters.SaddleMatrix = SaddleMatrix; 

IndexMatrix = load('MISE_IndexMatrix_Gauss4');
IndexMatrixNames = fieldnames(IndexMatrix); 
IndexMatrix = IndexMatrix.(IndexMatrixNames{1}); 
Parameters.IndexMatrix = IndexMatrix; 


%3. LC identifies the optimal parameter combination for a given goal,
%which is a function that evaluates the action of the dynamical system.
%That function must be defined for LC. 
%This function should be one which takes in the action for each transition
%in a dynamical system, and returns a single value. 
%For example, in the paper we use the occupancy of a single state for this
%value, as well as particular transition rates.  

LC_Functional = @(ActionArray, StableStatesInExistence)MISE_LCFunctional(ActionArray, StableStatesInExistence); 
Parameters.Functional = LC_Functional;  


%4. Stable State Identifier
%In the course of the simulation bifurcations can occur and two states can
%become one. A function is required to distinguish different stable 
%states. 
ID_Function = @(x)MISE_SSC(x, StableStates); 
Parameters.ID_Function = ID_Function; 


%5. Input the (controllable) Jacobian of the dynamical system. 
%This should take in the time, (t), state vector (x), and the set of
%parameters, and return the square matrix that is the Jacobian of the
%dynamical system. 

%Importantly, the Jacobian MUST be able to be vectorized - it must be able
%to be evaluated over a time series of different states of the system. 
%See the included example for details. 

%While this Jacobian is not technically required to run LC, the speed up 
%when it is present is immense (>100x), making it effectively required for 
%any realistic examples. 

Jacobian = @(t, x, P) MISE_Jacobian(t, x, P); 
Parameters.Jacobian = Jacobian;
Parameters.JacobianPresent = 1; %Set this to 0 if no analytic Jacobian is given. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LC Parameters
%Parameters to tune the main optimization step in LC.
%These parameters might need tuning for optimal performance. 
LC_Parameters = {}; 

LC_Parameters.MaxIter = 20; %Set the maximum number of iterations for LC. 
LC_Parameters.ObjectiveLimit = -1; %Set a tolerance for the objective function, defined above. 
LC_Parameters.ub = 1*ones(1, 4); %Bound parameter can be included to to maintain physiologically realistic parameter values. 
LC_Parameters.lb = 0*ones(1, 4); %lower bound for the parameters. 
%These bounds can be modified to achieve better convergence.


LC_Parameters.InitParams = [.5 .5 .5 .5]; %The initial parameter set for the model. 
%IMPORTANT: This parameter set should correspond to the stable state above.

Parameters.LC_Parameters = LC_Parameters;

%Make iterative figures? Note that this function is not modular and will
%need to be adapted to a different model, if desired. 
Parameters.MakeFigures = 0; 
 

%%%%%%%%%%%%%%%%%%%%%%%%%Bifurcation Constraint Function%%%%%%%%%%%%%%
%Bifurcation Constraint Function
%6.
%A nonlinear constraint should be given that gives conditions about how
%close each stable state can get to bifurcation. 

%Other constraints of both inequality and equality type, can be stipulated
%as well, in accordance with the documentation of fmincon. 


NonlinearConstraintParams.Function = Function;
NonlinearConstraintParams.Jacobian = Jacobian; 
NonlinearConstraintParams.ID_Function = ID_Function;
NonlinearConstraintParams.OriginalFP = StableStates;
NonlinearConstraintParams.OriginalParams = LC_Parameters.InitParams;
NonlinearConstraintParams.Tolerance = .35; %This value can be modified 
%to ensure proper functioning of the eigenvalue check procedure. 
 
ConstraintFunction = @(x)MISE_Constraint(x, NonlinearConstraintParams); 

%This mixing time constraint used in the example. 
%ConstraintFunction = @(x)MISE_Constraint_MixingTime(x, Parameters); 
Parameters.ConstraintFunction = ConstraintFunction; 
 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now execute LC. 
%With all of the parameters set up, the execution is easy. 


%EXECUTE!
OutParameters = OptimalBHControl(Parameters); 

%OutParameters is a struct which includes the following fields: 
OutParameters.ObjVal %What is the ultimate objective value achieved with
%LC. 
OutParameters.Params(:) %What is the set of parameter values that achieves
OutParameters.History.ObjectiveValues(:)
OutParameters.History.Params{:}%Give the full history of the optimization/



