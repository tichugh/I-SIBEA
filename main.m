function [Utility,Final_solution] = main(Parameters)

gen = Parameters(1).p;
Pop_Size = Parameters(2).p;
H = Parameters(3).p;
W_vector = Parameters(4).p;
Problem = Parameters(5).p; 
M = Parameters(6).p;
N_DI = Parameters(7).p; 



addpath(genpath('Public'));
addpath(genpath('Hypervolume'));
addpath(genpath('Decision_Making'));
addpath(genpath('Project_ASF'));
Real_PF = P_objective('true', Problem, M, 1000);
Ideal = min(Real_PF);
[Population,Boundary] = P_objective('init',Problem,M,Pop_Size);
V = size(Boundary,2);
FunctionValue = P_objective('value',Problem,M,Population);
FrontValue = P_sort(FunctionValue);
Loss = F_hypervolume(FunctionValue,FrontValue);
% Utility_vector=zeros(H+1,1);

N_DM = Pop_Size; % Maximum Number of solutions shown to the decision maker
HH=1;
gen_DM =0;

for i = 1:gen

    MatingPool = F_mating(Population,FrontValue,Loss);
    Offspring = P_generator(MatingPool,Boundary,'Real',Pop_Size);     
    Fitness = P_objective('value',Problem,M,Offspring);    
    
    Population = [Population;Offspring];
    FunctionValue = [FunctionValue;Fitness];
    [FrontValue,MaxFront] = P_sort(FunctionValue,'half');
    no_gen = i
    
    if i<=N_DI && HH ==1
        Loss = F_hypervolume(FunctionValue,FrontValue);
    else
        Loss = F_hypervolume_weight(FunctionValue,FrontValue,W_I,AA,RA,V);        
    end
    Next = zeros(1,Pop_Size);
    NoN = numel(FrontValue,FrontValue<MaxFront);
    Next(1:NoN) = find(FrontValue<MaxFront);
    Last = find(FrontValue==MaxFront);
    [~,Rank] = sort(Loss(Last),'descend');
    Next(NoN+1:Pop_Size) = Last(Rank(1:Pop_Size-NoN));
    Population = Population(Next,:);
    FrontValue = FrontValue(Next);
    Loss = Loss(Next);
    FunctionValue = FunctionValue(Next,:);
    if M==2
        scatter(FunctionValue(:,1),FunctionValue(:,2))
        axis([0,1.5,0,1.5]);
    elseif M==3
        scatter3(FunctionValue(:,1),FunctionValue(:,2),FunctionValue(:,3))
        axis([0,3.0,0,3.0,0,3.0]);
    else
        parallelcoords(FunctionValue);
        axis([1,M,0,2]);
    end
%     pause(0.1)

%% Final value using ASF
    if (i>=gen)
        display('Maximum generations are completed')
        if (H ~=0)
            r_rows = P_sort(FunctionValue,'first')==1;
            Population = Population(r_rows,:);
            FunctionValue = FunctionValue(r_rows,:);
            % project it to the Pareto front by solving ASF
            Auto = FunctionValue;
            Auto1 = bsxfun(@minus, Auto,Ideal);
            Auto1 = bsxfun(@times, Auto1,W_vector);
            Auto1 = max(Auto1,[],2);
            [~,Ind_AA] = min(Auto1);
            AA = [Population(Ind_AA,:),FunctionValue(Ind_AA,:)];
    %         [Utility_vector(HH,:),Final_solution] = project_ASF_new(AA,M,V,Boundary,W_vector,Problem,Ideal);
            [Utility,Final_solution] = project_ASF(AA,M,V,Boundary,W_vector,Problem,Ideal);
%             display('Maximum generations are completed')
        else
            Utility = [];
            Final_solution = [Population,FunctionValue];
           
        end
        if M==2
            scatter(Final_solution(:,V+1),Final_solution(:,V+2),'red')
            axis([0,1.5,0,1.5]);
        elseif M==3
            scatter3(Final_solution(:,V+1),Final_solution(:,V+2),Final_solution(:,V+3),'red')
            axis([0,3.0,0,3.0,0,3.0]);
        end
        legend('Final Solution(s)')
    end

    gen_DM =gen_DM + 1;
    %% Role of decision maker
    if gen_DM == N_DI && HH <= H
      
        display(['Role of decision maker-' 'Interaction Number -' num2str(HH)]);
%         [W_I,AA,RA,Utility]=Role_decision_maker([Population,FunctionValue],N_DM,M,V,W_vector,Ideal);
        [W_I,AA,RA]=Role_decision_maker([Population,FunctionValue],N_DM,M,V,W_vector,Ideal);
%         Utility_vector(HH,:) = Utility;
        HH = HH+1; 
        gen_DM = 0;
        N_DI = round(gen/(2*(H-1)));
    end

end 
end 









    











