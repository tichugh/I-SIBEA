function Loss = F_hypervolume(FunctionValue,FrontValue)
    [N,M] = size(FunctionValue);
    r = max(FunctionValue);
    Loss = zeros(1,N);
    Fronts = setdiff(unique(FrontValue),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontValue==Fronts(f));
        Fit = FunctionValue(Front,:);
        W = ones(1,size(Fit,1));
        TH = P_evaluate_hv_weight('nHV',Fit,r,W);
        temp = Fit;
        temp_W = W;
        
        for i = 1:size(temp,1)
            Fit(i,:)= [];
            W(i) = [];
            hv = P_evaluate_hv_weight('nHV',Fit,r,W);
            Loss(Front(i)) = TH - hv;
            Fit = temp;
            W = temp_W;
        end
        for i = 1 : M
            [~,Rank] = sortrows(FunctionValue(Front,i));
            Loss(Front(Rank(1))) = inf;
            Loss(Front(Rank(end))) = inf;
        end
    end
end

