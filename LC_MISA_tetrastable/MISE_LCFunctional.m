function [out, BHs, occupancies]=MISE_LCFunctional(BHs, FP)
%For our test functional, we optimize the occupancy of the third fixed
%point in the MISE system, assuming a noise level of epsilon = 0.02;
 
 
%Find those transitions which are faster when indirect. 
sz = size(BHs, 2);
% for i=1:sz
%     for j=1:sz
%         for k=1:sz
%             if actions(i,j) + actions(j,k) < actions(i,k)
%                 actions(i,k) = Inf; 
%             end
%         end
%     end
% end
 
BHs(logical(eye(sz))) = Inf; 
Epsilon = 0.35; 
rates=exp(-BHs/Epsilon);
% disp(rates)
rates(logical(eye(sz))) = -sum(rates,2); 
occupancies = null(rates')/sum(null(rates')); 
% occupancies = BHs/sum(BHs);
% disp(occupancies)
%This calculates the occupancies in a canonical way by finding the steady
%state of the continuous time Markov chain. 
 
%This is the occupancy of the second state. 
out = -sum(occupancies(FP==2)) ;
 
end
