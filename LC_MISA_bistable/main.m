%Define the parameter struct on which the algorithm will operate. 
Parameters = {}; 

Function  = @(t, x, P) MutualInhibitionSelfExcitation(t, x, P);
Parameters.Function = Function;

StableStates = load('MISE_StableStates2');
StableStateNames = fieldnames(StableStates); 
StableStates = StableStates.(StableStateNames{1}); 
Parameters.StableStates = StableStates; 

SaddleMatrix = load('MISE_SaddleMatrix_Gauss2');
SaddleMatrixNames = fieldnames(SaddleMatrix); 
SaddleMatrix = SaddleMatrix.(SaddleMatrixNames{1}); 
Parameters.SaddleMatrix = SaddleMatrix; 

IndexMatrix = load('MISE_IndexMatrix_Gauss2');
IndexMatrixNames = fieldnames(IndexMatrix); 
IndexMatrix = IndexMatrix.(IndexMatrixNames{1}); 
Parameters.IndexMatrix = IndexMatrix; 

LC_Functional = @(ActionArray, StableStatesInExistence)MISE_LCFunctional(ActionArray, StableStatesInExistence); 
Parameters.Functional = LC_Functional;  


ID_Function = @(x)MISE_SSC(x, StableStates); 
Parameters.ID_Function = ID_Function; 

Jacobian = @(t, x, P) MISE_Jacobian(t, x, P); 
Parameters.Jacobian = Jacobian;
Parameters.JacobianPresent = 1; %Set this to 0 if no analytic Jacobian is given. 

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

NonlinearConstraintParams.Function = Function;
NonlinearConstraintParams.Jacobian = Jacobian; 
NonlinearConstraintParams.ID_Function = ID_Function;
NonlinearConstraintParams.OriginalFP = StableStates;
NonlinearConstraintParams.OriginalParams = LC_Parameters.InitParams;
NonlinearConstraintParams.Tolerance = .35; 

ConstraintFunction = @(x)MISE_Constraint(x, NonlinearConstraintParams); 

%This mixing time constraint used in the example. 
%ConstraintFunction = @(x)MISE_Constraint_MixingTime(x, Parameters); 
Parameters.ConstraintFunction = ConstraintFunction; 
 
OutParameters = OptimalBHControl(Parameters); 

%OutParameters is a struct which includes the following fields: 
OutParameters.ObjVal %What is the ultimate objective value achieved with
%LC. 
OutParameters.Params(:) %What is the set of parameter values that achieves
OutParameters.History.ObjectiveValues(:)
OutParameters.History.Params{:}%Give the full history of the optimization/
