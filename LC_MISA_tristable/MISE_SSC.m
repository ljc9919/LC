function out = MISE_SSC(x, StableStates)
NumFP = size(StableStates, 2); 

[dist, out] = min(sum((repmat(x, [1, NumFP]) - StableStates).^2)); 
    
end
    