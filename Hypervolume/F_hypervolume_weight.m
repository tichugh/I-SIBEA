function Loss = F_hypervolume_weight(FunctionValue,FrontValue,W_I,AA,RA,V)

    [N,M] = size(FunctionValue);
    r = max(FunctionValue) + 0.1;
    Loss = zeros(1,N);
    Fronts = setdiff(unique(FrontValue),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontValue==Fronts(f));
        Fit = FunctionValue(Front,:);
        weight = ones(1,size(Fit,1));
        for i=1:size(Fit,1)
            my_temp = bsxfun(@le,Fit(i,:),AA(:,V+1:V+M));
            r_temp = find(all(my_temp==1,2));
            if (~isempty(r_temp))
                weight(i) = W_I(3,:);
            end
            my_temp2 = bsxfun(@gt,Fit(i,:),RA(:,V+1:V+M));
            r_temp2 = find(all(my_temp2==1,2));
            if (~isempty(r_temp2))
                weight(i) = W_I(1,:);
            end  
        end 
        TH = P_evaluate_hv_weight('nHV',Fit,r,weight);
        temp = Fit;
        temp2 = weight;
        for i = 1:size(temp,1)
            Fit(i,:)= [];
            weight(i) = [];
            hv = P_evaluate_hv_weight('nHV',Fit,r,weight);
            Loss(Front(i)) = TH - hv;
            Fit = temp;
            weight = temp2;
        end

    end
end

