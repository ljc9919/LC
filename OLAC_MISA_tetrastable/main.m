%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OptimalLeastActionControl_Top.m


%MATLAB implementation of the computational method Optimal Least Action
%Control to accompany the paper "Control of Stochastic and Induced
%Switching in Biophysical Complex Networks". 
%Daniel K. Wells, William L. Kath & Adilson E. Motter
%Northwestern University 
%Copyright ? 2015 All rights reserved.
%Ver 1.0
%Contact Information: Daniel K. Wells: dannykwells@gmail.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%License information for this package (and those packages/scripts above, including: 
%OptimalLeastActionControl_Top.m , OptimalLeastActionControl.m,
%MISE_Constraint.m,%MISE_Constraint_MixingTime.m, MISE_Jacobian.m,
%MISE_OLACFunctional.m, MISE_SSC.m, MutualInhibitionSelfExcitation.m):


%Redistribution and use in source and binary forms, with or without
%modification, are permitted provided that the following conditions are
%met: 
%*Redistributions of source code must retain the the license information,
%this list of conditions, and the following disclaimer in the
%documentation and/or materials provided with distribution. 
%*Redistribution in binary form must reproduce the above license
%information, this list of conditions, and the following disclaimer in the
%documentation and/or materials provided with distribution. 

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%TOP LEVEL: DEFINE FUNCTIONS AND PARAMETERS. 
%For basic model exploration and least action control. 
%Advanced users can modify functions/parameters in the other scripts. 
%This code will allow one to use OLAC in a standard way. Complicated
%objective function or constraints, very large systems, or other
%complications might require modifications. 


%This setup will work "out of the box" on the included functions. 
%Changes will need to be made to have it work on a different
%problem. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Define the parameter struct on which the algorithm will operate. 
Parameters = {}; 



%1. Define the (controllable dynamical system). 
%This should be a function that takes in three parameters, in this order: 
%(t, x, P), where t is the time (for non-autonomous dynamical systems), x
%is the state vector, and P is the parameter vector for the system. 


Function  = @(t, x, P) MutualInhibitionSelfExcitation(t, x, P);
Parameters.Function = Function;



%2. The states of the system must be previously identified to run OLAC. These
%should be given as an n-row by m-column array, where n is the number of
%dimensions of the dynamical system under consideration, and m is the
%number of stable states. 
%should be inputted here:

StableStates = load('MISE_StableStates4');
StableStateNames = fieldnames(StableStates); 
StableStates = StableStates.(StableStateNames{1}); 
Parameters.StableStates = StableStates; 


%3. OLAC identifies the optimal parameter combination for a given goal,
%which is a function that evaluates the action of the dynamical system.
%That function must be defined for OLAC. 
%This function should be one which takes in the action for each transition
%in a dynamical system, and returns a single value. 
%For example, in the paper we use the occupancy of a single state for this
%value, as well as particular transition rates.  

OLAC_Functional = @(ActionArray, StableStatesInExistence)MISE_OLACFunctional(ActionArray, StableStatesInExistence); 
Parameters.Functional = OLAC_Functional;  


%4. Stable State Identifier
%In the course of the simulation bifurcations can occur and two states can
%become one. A function is required to distinguish different stable 
%states. 40
ID_Function = @(x)MISE_SSC(x, StableStates); 
Parameters.ID_Function = ID_Function; 


%5. Input the (controllable) Jacobian of the dynamical system. 
%This should take in the time, (t), state vector (x), and the set of
%parameters, and return the square matrix that is the Jacobian of the
%dynamical system. 

%Importantly, the Jacobian MUST be able to be vectorized - it must be able
%to be evaluated over a time series of different states of the system. 
%See the included example for details. 

%While this Jacobian is not technically required to run OLAC, the speed up 
%when it is present is immense (>100x), making it effectively required for 
%any realistic examples. 

Jacobian = @(t, x, P) MISE_Jacobian(t, x, P); 
Parameters.Jacobian = Jacobian;
Parameters.JacobianPresent = 1; %Set this to 0 if no analytic Jacobian is given. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OLAC Parameters
%Parameters to tune the main optimization step in OLAC.
%These parameters might need tuning for optimal performance. 
OLAC_Parameters = {}; 

OLAC_Parameters.MaxIter = 20; %Set the maximum number of iterations for Olac. 
OLAC_Parameters.ObjectiveLimit = -1; %Set a tolerance for the objective function, defined above. 
OLAC_Parameters.ub = 1*ones(1, 4); %Bound parameter can be included to to maintain physiologically realistic parameter values. 
OLAC_Parameters.lb=0*ones(1, 4); %lower bound for the parameters. 
%These bounds can be modified to achieve better convergence.


OLAC_Parameters.InitParams =[.4 .6 .4 .6]; %The initial parameter set for the model. 
%IMPORTANT: This parameter set should correspond to the stable state above.

Parameters.OLAC_Parameters = OLAC_Parameters;

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
NonlinearConstraintParams.OriginalFP = StableStates;
NonlinearConstraintParams.OriginalParams = OLAC_Parameters.InitParams;
NonlinearConstraintParams.Tolerance = .35; %This value can be modified 
%to ensure proper functioning of the eigenvalue check procedure. 
 
ConstraintFunction = @(x)MISE_Constraint(x, NonlinearConstraintParams); 

%This mixing time constraint used in the example. 
%ConstraintFunction = @(x)MISE_Constraint_MixingTime(x, Parameters); 
Parameters.ConstraintFunction = ConstraintFunction; 
 







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Now execute OLAC. 
%With all of the parameters set up, the execution is easy. 


%EXECUTE!
OutParameters = OptimalLeastActionControl(Parameters); 

%OutParameters is a struct which includes the following fields: 
OutParameters.ObjVal %What is the ultimate objective value achieved with
%OLAC. 
OutParameters.Params(:) %What is the set of parameter values that achieves
OutParameters.History.ObjectiveValues(:)
OutParameters.History.Params{:}%Give the full history of the optimization/
