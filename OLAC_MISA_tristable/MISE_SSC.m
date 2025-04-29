%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MISE_SSC.m
%A function to classify new stable states of the system according to which 
%stable state that point originated from, i.e., is a continuation of. 
%We use a shortest distance metric in this case. 
%Other more biologically motivated examples can be used as well. 

%Author: Daniel K. Wells (©) 2015
%Ver 1.0
%Email: dannykwells@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function out = MISE_SSC(x, StableStates)
NumFP = size(StableStates, 2); 

[dist, out] = min(sum((repmat(x, [1, NumFP]) - StableStates).^2)); 
    
end
    