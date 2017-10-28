function [Utility,my_solution] = project_ASF(Solutions,M,V,Boundary,W,Problem,Ideal)
lb = Boundary(2,:);
ub = Boundary(1,:);

Utility = zeros(size(Solutions,1),1);
my_solution = zeros(size(Solutions,1),V+M);
for i=1:size(Solutions,1)
    x0(1,1) = 0;
    x0(1,2:V+1)=Solutions(i,1:V);
    x0=x0';
    zref = Solutions(i,V+1:V+M);
    options = optimoptions(@fmincon,'Algorithm','sqp');
    x=fmincon(@(x)myfun(x),x0,[],[],[],[],lb,ub,@(x)constraints(x,M,zref,W,Problem,V),options);
    x=x';
    fun = P_objective('value',Problem,M,x(2:V+1));
    my_solution(i,:) = [x(2:end),fun];
    Auto = my_solution(i,V+1:V+M);
    Auto1 = bsxfun(@minus, Auto,Ideal);
    Auto1 = bsxfun(@times, Auto1,W);
    Auto1 = max(Auto1,[],2);
    Utility(i,:) = min(Auto1);
end
end

function asf=myfun(x)
    asf = x(1);
end

function [c,ceq] = constraints(x,M,zref,W,Problem,V)
x = x';
f = P_objective('value',Problem,M,x(2:V+1));
% f = obj_fun(x(2:V+1),M,V);   
c = zeros(1,M);
for k = 1:M
    c(k) =  W(k)*(f(k)-zref(k)) - x(1);
end
ceq = [];
end