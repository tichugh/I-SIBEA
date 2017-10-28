function [W_I,AA,RA,Utility] = Role_decision_maker(parents,N_DM,M,V,weight_vector,Ideal)
%% Here preferred and non-preferred solutions are selected based on the automated version using Chebyshve function values.
% To incoportate the real decision maker, one needs to modify this file
% as per the article

Rrows = P_sort(parents(:,V+1:V+M),'first')==1;
non_dominated_parents = parents(Rrows,:);
r = max(non_dominated_parents(:,V+1:V+M))+0.0001;
% ask the decision maker to partition the solutions into AA and RA
[~, ic]= unique(non_dominated_parents(:,V+1:V+M),'rows');  
non_dominated_parents = non_dominated_parents(ic,:);


    
%%    Automated version for Indices_AA
Auto = non_dominated_parents(:,V+1:V+M);
Auto1 = bsxfun(@minus, Auto,Ideal);
Auto1 = bsxfun(@times, Auto1,weight_vector);
Auto1 = max(Auto1,[],2);

[Utility,Indices_AA] = min(Auto1);
AA = non_dominated_parents(Indices_AA,:);
non_dominated_parents(Indices_AA,:) = [];
RA = non_dominated_parents;

W1 = ones(1,size(AA,1));
W2 = ones(1,size(RA,1));
hypervolume_Pr= P_evaluate_hv_weight('nHV',AA(:,V+1:V+M),r,W1);
hypervolume_Do= P_evaluate_hv_weight('nHV',RA(:,V+1:V+M),r,W2);
    
I_D = 0;
I_I = 1;
I_P = 1 + (hypervolume_Do/hypervolume_Pr); % Possibility 1
W_I = [I_D I_I I_P]';
   
