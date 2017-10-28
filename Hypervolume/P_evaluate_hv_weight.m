function Metric = P_evaluate_hv_weight(Operation,FunctionValue,RefPoint,W,varargin)


    [N,M] = size(FunctionValue);
    
    
    switch Operation
        case 'HV'
            
            FunctionValue(sum(FunctionValue-repmat(RefPoint,N,1)<=0,2)<M,:) = [];
            N = size(FunctionValue,1);
            if isempty(FunctionValue)
                Metric = 0;
                return;
            elseif M < 4
                pl     = sortrows(FunctionValue,1);
                S(1,1) = {1};
                S(1,2) = {pl};
                for k = 1 : M-1
                    S_ = {};
                    j  = size(S,1);
                    for l = 1 : j
                        S_l    = cell2mat(S(l,2));
                        S_temp = Slice(S_l,k,RefPoint);
                        p      = size(S_temp,1);
                        for q = 1 : p
                            cell_(1,1) = {cell2mat(S_temp(q,1))*cell2mat(S(l,1))};
                            cell_(1,2) = S_temp(q,2);
                            S_         = Add(cell_,S_);
                        end
                    end
                    S = S_;
                end
                Vol = 0;
                for l = 1 : size(S,1)
                    p   = Head(cell2mat(S(l,2)));
                    num = abs(p(M)-RefPoint(M));
                    Vol = Vol+cell2mat(S(l,1))*num*W(l);
                end
            else
%                 if nargin > 4
                    k_samples = varargin{1};
                    k_samples = floor(k_samples);
                    k_T        = varargin{2};        
%                 else
%                     k_T = 1000;
%                     k_samples = k_T;
%                 end
                
                MaxValue = RefPoint;
                MinValue = min(FunctionValue,[],1);
%                 Samples  = repmat(MinValue,k,1)+rand(k,M).*repmat((MaxValue-MinValue),k,1);
                
                Samples  = repmat(MinValue,k_samples,1)+rand(k_samples,M).*repmat((MaxValue-MinValue),k_samples,1);
                Domi     = false(1,k_samples);
                for i = 1 : N
                    Domi(sum(repmat(FunctionValue(i,:),k_samples,1)-Samples<=0,2)==M) = true;
%                     
                end
                Vol = prod(MaxValue-MinValue)*sum(Domi)/k_T;
            end
            Metric = Vol;           
        case 'nHV'
            if M>=4
                if nargin > 4
                    k_samples = varargin{1};
                    k_T        = varargin{2};  
                else
                    k_T = 10000;
                    k_samples = k_T;
                end
            end
            Metric = P_evaluate_hv_weight('HV',FunctionValue,RefPoint,W,k_samples,k_T)/abs(prod(max(RefPoint,[],1)));
    end
end

function S = Slice(pl,k,refPoint)
    p  = Head(pl);
    pl = Tail(pl);
    ql = [];
    S  = {};
    while ~isempty(pl)
        ql  = Insert(p,k+1,ql);
        p_  = Head(pl);
        cell_(1,1) = {abs(p(k)-p_(k))};
        cell_(1,2) = {ql};
        S  = Add(cell_,S);
        p  = p_;
        pl = Tail(pl);
    end
    ql = Insert(p,k+1,ql);
    cell_(1,1) = {abs(p(k)-refPoint(k))};
    cell_(1,2) = {ql};
    S  = Add(cell_,S);
end

function ql = Insert(p,k,pl)
    flag1 = 0;
    flag2 = 0;
    ql    = [];
    hp    = Head(pl);
    while ~isempty(pl) && hp(k) < p(k)
        ql = Append(hp,ql);
        pl = Tail(pl);
        hp = Head(pl);
    end
    ql = Append( p , ql );
    m  = length(p);
    while ~isempty(pl)
        q = Head(pl);
        for i = k : m
            if p(i) < q(i)
                flag1 = 1;
            else
                if p(i) > q(i)
                    flag2 = 1;
                end
            end
        end
        if ~(flag1 == 1 && flag2 == 0)
            ql = Append(Head(pl),ql);
        end
        pl = Tail(pl);
    end  
end

function p = Head(pl)
    if isempty(pl)
        p = [];
    else
        p = pl(1,:);
    end
end

function ql = Tail(pl)
    m = size(pl,1);
    if m == 0 || m == 1
        ql = [];
    else
        ql = pl(2:m,:);
    end
end

function S_ = Add(cell_,S)
    n = size(S,1);
    m = 0;
    for k = 1 : n
        if isequal(cell_(1,2),S(k,2))
            S(k,1) = {cell2mat(S(k,1))+cell2mat(cell_(1,1))};
            m = 1;
            break;
        end
    end
    if m == 0
        S(n+1,:) = cell_(1,:);
    end
    S_ = S;     
end

function pl = Append(p,ql)
    m  = size(ql,1);
    ql(m+1,:) = p;
    pl = ql;
end