% This function implements the Greedy match routine to match the leftover
% supplies and demands from the GTTransport
% implementations. 

function greedyCapacity = greedyMatch(n, gsupplies, gdemands, CostGreedy)

% Greedy matching
greedyCapacity = zeros(n);

% Sort using quick sort
[~,I] = sort(CostGreedy, 2);

for i = 1:n
    j = 1;
    while gsupplies(i) > 0 && j <= n
        in2 = I(i,j);
        units = min(gsupplies(i),gdemands(in2));
        greedyCapacity(i,in2) = units;
        gsupplies(i) = gsupplies(i) - units;
        gdemands(in2) = gdemands(in2) - units;
        CostGreedy(i,in2) = Inf;
        j = j+1;
    end
end

end

