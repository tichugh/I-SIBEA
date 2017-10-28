%% This is the code for I-SIBEA implemented by Tinkle Chugh
% Please cite the following paper if you use the code:
% Chugh et al. An interactive simple indicator based evolutionary algorithm (I-SIBEA) 
% for multiobjective optimization problems. In: 8th Int. Conf. on Evolutionary Multi-Criterion Optimization (EMO'2015), 
% Springer International Publishing, 277-291, 2015 
% Contact: tinkle.chugh@gmail.com
%%
function main_I_SIBEA
%% Here preferred and non-preferred solutions are selected based on the automated version using Chebyshve function values.
% To incoportate the real decision maker, one needs to modify this file
% "Role_decision_maker" as per the article

% Recommended to use the small population size, this version is implemented
% upto five objectives

%% The code can be used as apriori if number of interactions (H)=1  
% and as post-eriori if the number of interactions (H)=0

% mode = 'a-priori';
% mode = 'interactive';
mode = 'post-eriori';
w_vector = [1,1,1,1]; %% weights to calculate the scalarizing function value
Problem = 'DTLZ2';
M = 4; % number of objectives
%% Parameters 
Generations=200;
Pop_Size=10;
MaxRun = 1;

switch mode
    case {'a-priori'}
        No_Interactions = 1; % Number of interactions
        N_DI = round(Generations/3); % number of generations before the first interaction
    case {'interactive'}
        No_Interactions = 4; % Number of interactions
        N_DI = round(Generations/3); % number of generations before the first interaction
    case {'post-eriori'}
        No_Interactions = 0; % Number of interactions
        N_DI = Generations; % number of generations before the first interaction
end



%% Main loop
for i=1:length(No_Interactions)
    H = No_Interactions(i);    
    for run_no=1:MaxRun
        
        Parameters = struct('p',{Generations,Pop_Size,H,w_vector,Problem, M,N_DI});
        [Utility_value,Final_solution] = main(Parameters)  
    end
end