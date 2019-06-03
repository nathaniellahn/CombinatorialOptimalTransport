% This function implements the main functionality for the 1-feasibility transportation
% algorithm.

% Gabow Tarjan implementation for Optimal Transport Problem
function [iterationCount, APLengths, capacity] =  GTTransport(C, supplies, demands, n)

    APLengths = 0;%The total length of all augmenting paths computed thus far.
    % Source and sink nodes
    source = n*2+1;
    sink = n*2+2;

    % Set BFree and set AFree of vertices
    BFree = ones(1,n);
    AFree = ones(1,n);

    % Appending sink and source values to the cost matrix
    row = zeros(1,n*2);
    C = [C;row;row];
    col = zeros(1,n*2+2)';
    C = [C col];
    C = [C col];

    % Dual weights matrix
    dualWeights = zeros(1,(n*2)+2);

    % Deficiencies
    deficiency(1:n) = supplies;
    deficiency(n+1:n*2) = demands;

    for i=1:n
        if ~deficiency(i)
            BFree(i)=0;
        end
    end

    for i=n+1:n*2
        if ~deficiency(i)
            AFree(i-n)=0;
        end
    end

    % Edge Capacities
    capacity = zeros(n*2+2);
    for i = n+1:n*2
        for j = 1:n
            capacity(j,i) = min(deficiency(j),deficiency(i));
        end
    end

    % Initializing count for iteration and the cardinality of the matching
    iterationCount = 0;

    % One scale of the Gabow Tarjan algorithm
    % Continue while the matching is not perfect
    while any(BFree(:))

        % Incrementing the iteration count
        iterationCount = iterationCount + 1;

        % Creating residual graph
        residualGraph = zeros(n*2+2);
        for i=1:n
            for j=n+1:n*2
                if capacity(j,i)>0
                    residualGraph(j,i) = 1;
                end
            end
            if BFree(i)
                residualGraph(source,i) = 1;
            end
            if AFree(i)
                residualGraph(i+n,sink) = 1;
            end
        end

        for i = n+1:n*2
            for j = 1:n
                if capacity(j,i)>0
                    residualGraph(j,i) = 1;
                end
            end
        end

        % Dijkstra's algorithm to perform the Hungarian search
        lv = Inf(1, n*2+2);
        lv(source) = 0;
        shortestpathset = false(1, n*2+2);

        for i = 1:n*2+2
            minIndex = DijkstraMinDistance(lv,shortestpathset, n);
            shortestpathset(minIndex) = true;
            for j = 1:n*2+2
                if ~shortestpathset(j) && residualGraph(minIndex,j)
                    if j<n+1
                        slack = abs(C(minIndex,j)-dualWeights(minIndex)-dualWeights(j));
                    else
                        slack = abs(C(minIndex,j)+1-dualWeights(minIndex)-dualWeights(j));
                    end
                    if minIndex==source
                        slack = 0;
                    end
                    if j==sink
                        slack = 0;
                    end
                    if lv(j) > lv(minIndex) + slack
                        lv(j) = lv(minIndex) + slack;
                    end
                end
            end
            if minIndex > n && minIndex ~= source && minIndex~= sink && AFree(minIndex-n)
                break;
            end
        end

        % Update dual weights
        for i=1:n*2+2
            if lv(i)<lv(sink)
                if i<n+1
                    dualWeights(i) = dualWeights(i) + lv(sink) - lv(i);
                elseif i>n && i<n*2+1
                    dualWeights(i) = dualWeights(i) - lv(sink) + lv(i);
                end
            end
        end

        % Create Admissible graph
        admissibleGraph = zeros(n*2);
        for i=1:n
            for j=n+1:n*2
                if residualGraph(j,i) && (abs(C(i,j) - dualWeights(i) - dualWeights(j)) == 0)
                    admissibleGraph(j,i) = 1;
                end
            end
        end
        for i = n+1:n*2
            for j = 1:n
                if residualGraph(j,i) && ((C(i,j) + 1 - dualWeights(i) - dualWeights(j)) == 0)
                    admissibleGraph(j,i) = 1;
                end
            end
        end

        % Partial DFS to get the maximal set of vertex-disjoint augmenting paths
        % Find augmenting paths for all the free vertices of set B
        augmentingPathVertices = zeros(1,n*2);
        for vertex=1:n
            vertexVisited = zeros(1,n*2);
            vertexVisited(1:n) = n;
            if BFree(vertex)

                % Keep finding augmenting paths for this vertex till the capacity is 0
                % or no more paths can be find
                while deficiency(vertex) > 0 && vertexVisited(vertex) < n*2


                    AugPathVerIndex = 1;
                    [augmentingPathVertices, vertexVisited, admissibleGraph, AugPathVerIndex] = DFSUtil(vertex, augmentingPathVertices, n, admissibleGraph, AFree, AugPathVerIndex, vertexVisited, capacity);

                    % Add the path to the maximal set of augmenting paths
                    [~,lengthOfPath] = size(augmentingPathVertices);
                    APLengths = APLengths + lengthOfPath;

                    if AugPathVerIndex > 1
                        % Calculate the initial flow
                        val = min(deficiency(augmentingPathVertices(1)), deficiency(augmentingPathVertices(AugPathVerIndex)));

                        % Calculate the flow and augment along it
                        for j = 1:AugPathVerIndex-1
                            vertex1 = augmentingPathVertices(j);
                            vertex2 = augmentingPathVertices(j+1);
                            val = min(val, capacity(vertex1,vertex2));
                        end

                        % Augment along the path
                        for j = 1:AugPathVerIndex-1
                            vertex1 = augmentingPathVertices(j);
                            vertex2 = augmentingPathVertices(j+1);
                            capacity(vertex1,vertex2) = capacity(vertex1,vertex2) - val;
                            capacity(vertex2,vertex1) = capacity(vertex2,vertex1) + val;
                            if capacity(vertex1,vertex2) > 0
                                vertexVisited(vertex1) = vertex2 -1;
                            end
                        end

                        % Change the deficiency at start node
                        deficiency(vertex) = deficiency(vertex) - val;
                        if deficiency(vertex) == 0
                            BFree(vertex) = 0;
                        end

                        % Change the deficiency at end node
                        deficiency(augmentingPathVertices(AugPathVerIndex)) = deficiency(augmentingPathVertices(AugPathVerIndex)) - val;
                        if deficiency(augmentingPathVertices(AugPathVerIndex)) == 0
                            AFree(augmentingPathVertices(AugPathVerIndex)-n) = 0;
                        end

                    end
                end
            end
        end

    end
end

% Utility function for Dijkstra's algorithm
function minIndex = DijkstraMinDistance(lv, shortestpathset, n)
minValue = Inf;
minIndex = 1;
for v = 1:n*2+2
    if ~shortestpathset(v) && lv(v) < minValue
        minValue = lv(v);
        minIndex = v;
    end
end
end

function [augmentingPathVertices, vertexVisited, admissibleGraph, AugPathVerIndex] = DFSUtil(vertex, augmentingPathVertices, n, admissibleGraph, AFree, AugPathVerIndex, vertexVisited, capacity)

    stack = vertex;
    while stack ~= 0
        vertex = stack;
        stack = 0;
        range_var1 = vertexVisited(vertex)+1;
        if vertex<n+1
            range_var2 = n*2;
        else
            range_var2 = n;
        end
        for i = range_var1:range_var2
            vertexVisited(vertex) = i;
            if admissibleGraph(vertex,i)
                if capacity(vertex,i) > 0
                    stack = i;
                    augmentingPathVertices(AugPathVerIndex) = vertex;
                    AugPathVerIndex = AugPathVerIndex + 1;
                    break;
                else
                    admissibleGraph(vertex,i) = 0;
                end
            end
        end
        if vertex>n && AFree(vertex-n)
            if augmentingPathVertices(AugPathVerIndex-1)~=vertex
                augmentingPathVertices(AugPathVerIndex) = vertex;
            end
            return;
        else
            if stack ~= 0
                continue;
            else
                if AugPathVerIndex>1
                    AugPathVerIndex = AugPathVerIndex - 1;
                    backtrack = augmentingPathVertices(AugPathVerIndex);
                    admissibleGraph(backtrack,vertex) = 0;
                    stack = backtrack;
                    augmentingPathVertices = augmentingPathVertices(augmentingPathVertices~=backtrack);
                end
            end
        end
    end
end












