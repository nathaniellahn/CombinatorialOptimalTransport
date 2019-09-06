% Gabow Tarjan implementation for integral Optimal Transport Problem.
% Input format: 
% ~ supplies is a [1 x n] column vector of integers that describes the 
% supplies at each vertex of type "B".
% ~ demands is a [1 x n] column vector of integers that describes the 
% demands at each vertex of type "A".
% This code requires that the total supplies are at most the total demands.
% This ensures that all supplies are matched by the end of the algorithm.
% ~ C : An [n x n] matrix of costs.
% ~ n : The problem size; specifically, the number of supply and demand 
% vertices (which are required to be the same).

function [iterationCount, APLengths, capacity] =  GTTransport(C, supplies, demands, n)
    %Check that mapping correctly ensures that supplies are less than
    % demands.
    assert(sum(supplies) <= sum(demands)); 
    
    %Cost matrix from B to A. 
    %Vertices of A are rows. Vertices of B are columns
    CBA = C';  
    
    %Cost matrix from A to B. 
    %Vertices of B are rows. Vertices of A are columns
    %While CBA could also be used, this allows column-order traversal to 
    %minimize cache misses when exploring neighbors
    CAB = C;
    
    APLengths = 0;%The total length of all augmenting paths computed thus far.
    % Source and sink nodes

    % Set BFree and set AFree of vertices
    BFree = ones(1,n);
    AFree = ones(1,n);
    
    % Dual weights matrix
    dualWeights = zeros(1,2*n);

    % Deficiencies
    deficiencyB = supplies;
    deficiencyA = demands;

    for i=1:n
        if ~deficiencyB(i)
            BFree(i)=0;
        end
        if ~deficiencyA(i)
            AFree(i) = 0;
        end
    end

    % Edge Capacities
    
    % Like all the other 2D matrices used, edges from A -> B have A as
    % columns.
    % Likewise, B -> A have vertices of B as columns.
    
    capacityAB = zeros(n);
    capacityBA = zeros(n);
    for j = 1:n %For each column
        for i = 1:n %For each row
            capacityBA(i,j) = min(deficiencyB(j),deficiencyA(i));
        end
    end

    % Initializing count for iteration and the cardinality of the matching
    iterationCount = 0;

    %Preallocate several arrays
    % As with other 2D matrices, the tail of the edge is the column.
    % The head is the row.
    slackAB = zeros(n);
    slackBA = zeros(n);
    vertexVisited = zeros(1,2*n);
    augmentingPathVertices = zeros(1,n*2);
    
    % Compute initial slacks
    for i=1:n %For each column
        for j=1:n  %For each row
            % First, assume i in A and j in B
            % This is an edge directed from A to B,
            % so the slack is given by the below expression
            slackAB(j,i) = dualWeights(i + n) + dualWeights(j) - CAB(j,i);

            % Next, assume i in B and j in A.
            % This is an edge from B to A.
            slackBA(j,i) = CBA(j,i) + 1 - dualWeights(i) - dualWeights(j + n);
        end
    end
    
    %The main phase loop.
    %Continue while the matching is not perfect
    %Each iteration runs in O(n^2) time 
    %(assuming augmenting path lengths are sufficienctly small, which
    %seems to be a good assumption).
    while any(BFree(:))
        % Incrementing the iteration count
        iterationCount = iterationCount + 1;
        
        % Dijkstra's algorithm to perform the Hungarian search
        %Distances to vertices measured in total slack from nearest
        %Free vertex of B.
        lv = Inf(1, n*2);
        
        %Distance to every free vertex of B is 0
        for i = 1:n
            if BFree(i)
                lv(i) = 0;
            end
        end
        
        %The set of vertices that have been reached in the Dijkstra search
        shortestpathset = false(1, n*2);
        
        % Record the minimum distance to a free vertex. 
        % -1 means no AP found.
        distAF = -1;
        
        %Main Dijkstra loop. Each iteration adds a vertex to the 
        %shortest path tree.
        %Will break out early if a free vertex of A is found.
        for i = 1:n*2
            %The next vertex to add to the shortest path tree.
            %Such distances are final, and will not be changed.
            % If -1 is returned, then all vertices are in shortestpathset
            % and the search should be terminated
            % however, this should never happen due to the for loop
            % condition.
            minIndex = DijkstraMinDistance(lv,shortestpathset, n);
            
            %Vertex added to shortest path tree
            shortestpathset(minIndex) = true; 
            
            if minIndex <= n
                %A vertex b of type B was added to the shortest path tree.
                %b = minIndex.
                %For each neighbor of b (iterate through column of BA
                %matrix), update distance to neighbors
                for a = 1:n 
                    if capacityBA(a,minIndex) > 0
                        slack = slackBA(a, minIndex);
                        if lv(a + n) > lv(minIndex) + slack
                            lv(a + n) = lv(minIndex) + slack;
                        end
                    end
                end
            elseif minIndex > n
                %A vertex of type A was added to the shortest path tree.
                a = minIndex - n;
                if AFree(a)
                    distAF = lv(minIndex);
                    break; % Dijkstra found an augmenting path
                end
                %For each neighbor of a (iterate through column of AB
                %matrix), update distance to neighbors
                for b = 1:n
                    if capacityAB(b,a) > 0
                        slack = slackAB(b, a); 
                        if lv(b) > lv(minIndex) + slack
                            lv(b) = lv(minIndex) + slack;
                        end
                    end
                end    
            end
        end
                
        %Since there are fewer supply vertices than demand vertices,
        %and the main loop only runs if a free supply vertex remains,
        %there should be no reason an augmenting path is not found
        %assuming all is set up properly.
         
        % Update dual weights
        for i=1:n*2
            if lv(i) < distAF
                if i <= n  %i is a vertex of B. Increase dual weight
                    dualWeights(i) = dualWeights(i) + distAF - lv(i);
                else % i > n, and is a vertex of A. Decrease dual weight
                    dualWeights(i) = dualWeights(i) - distAF + lv(i);
                end
            end
        end

        % Update slacks after dual adjustment, and compute admissible graph
        % for the upcoming DFS portion
        
        for i=1:n %For each column
            for j=1:n  %For each row
                % First, assume i in A and j in B
                % This is an edge directed from A to B,
                % so the slack is given by the below expression
                
                slack = dualWeights(i + n) + dualWeights(j) - CAB(j,i);
                slackAB(j,i) = slack;
                
                % Next, assume i in B and j in A.
                % This is an edge from B to A.          
                slack = CBA(j,i) + 1 - dualWeights(i) - dualWeights(j + n);
                slackBA(j,i) = slack;
            end
        end
        
        % Next, conduct several partial DFS searches to get the maximal set
        % of vertex-disjoint admissible augmenting paths.
        % A DFS search is initiated separately from each free vertex of B
        
        %This is used by the DFS procedure to track the
        %largest explored neighbor index of every vertex.
        %Following our convention, 1:n -> B and n+1:2*n -> A.
        %n is one less than the index of the first vertex of A
        vertexVisited(1:n) = n; 
        %0 is one less than the first vertex of B.
        vertexVisited(n + 1: 2*n) = 0;
        
        for vertex=1:n
            %For each free vertex
            if BFree(vertex)
                % Keep finding augmenting paths from this vertex
                % until either the vertex is no longer free
                % or no more paths can be found
                while deficiencyB(vertex) > 0 && vertexVisited(vertex) < n*2
                    %Initiate a partial DFS from this free vertex.
                    %~augmentingPathVertices is an array of vertices
                    %that describe an augmenting path.
                    %Note that everything past AugPathVerIndex is garbage.
                    %~AugPathVerIndex is the index of the last vertex of
                    %the augmenting path.
                    [augmentingPathVertices, vertexVisited, AugPathVerIndex] = DFSUtil(vertex, augmentingPathVertices, n, slackAB, slackBA, AFree, vertexVisited, capacityAB, capacityBA);
                    
                    if AugPathVerIndex <= 0
                        %No augmenting path found.
                        break;
                    else %AugPathVerIndex > 0                        
                        % Add the path to the maximal set of augmenting paths
                        % and update all relevant values.
                        
                        %Update total AP lengths for record keeping.
                        %Not used in actual algorithm.
                        APLengths = APLengths + AugPathVerIndex - 1;
                        
                        % Calculate the maximum flow the augmenting path
                        % can carry. This is the minimum of:
                        % 1.) The unmatched supplies at the starting free
                        % vertex of B.
                        % 2.) The unmatched demands at the ending free
                        % vertex of A
                        % 3.) The smallest capacity edge on the augmenting
                        % path itself.
                        % We compute and store this value in "val" below.
                        val = min(deficiencyB(augmentingPathVertices(1)), deficiencyA(augmentingPathVertices(AugPathVerIndex) - n));
                        for j = 1:AugPathVerIndex-1
                            vertex1 = augmentingPathVertices(j);
                            vertex2 = augmentingPathVertices(j+1);
                            if vertex1 > n
                                %Edge is directed from A to B
                                val = min(val, capacityAB(vertex2, vertex1 - n));
                            else
                                %Edge is directed from B to A
                                val = min(val, capacityBA(vertex2 - n, vertex1));
                            end
                        end

                        % Augment along the path
                        for j = 1:AugPathVerIndex-1
                            vertex1 = augmentingPathVertices(j);
                            vertex2 = augmentingPathVertices(j+1);
                            if vertex1 > n
                                %Edge is directed from A to B
                                capacityAB(vertex2, vertex1 - n) = capacityAB(vertex2, vertex1 - n) - val;
                                capacityBA(vertex1 - n, vertex2) = capacityBA(vertex1 - n, vertex2) + val;
                                if capacityAB(vertex2, vertex1 - n) > 0
                                    %Allow edge to be reused on future
                                    %augmenting paths.
                                    vertexVisited(vertex1) = vertex2 - 1; 
                                end
                            else
                                %Edge is directed from B to A
                                capacityBA(vertex2 - n, vertex1) = capacityBA(vertex2 - n, vertex1) - val;
                                capacityAB(vertex1, vertex2 - n) = capacityAB(vertex1, vertex2 - n) + val;
                                
                                if capacityBA(vertex2 - n, vertex1) > 0
                                    %Allow edge to be reused on future
                                    %augmenting paths.
                                    vertexVisited(vertex1) = vertex2 - 1; 
                                end
                            end
                        end %End augmentation for loop.
                        
                        % Next, we need to update the deficiencies
                        % of the endpoints of the path.
                        % Note that all other vertices' deficiencies
                        % are unchanged.
                        
                        % Change the deficiency at start node
                        deficiencyB(vertex) = deficiencyB(vertex) - val;
                        if deficiencyB(vertex) == 0
                            BFree(vertex) = 0;
                        end

                        % Change the deficiency at last node of path
                        last = augmentingPathVertices(AugPathVerIndex) - n;
                        deficiencyA(last) = deficiencyA(last) - val;
                        if deficiencyA(last) == 0
                            AFree(last) = 0;
                        end
                    end %End if/else checking whether an augmenting path
                        %was found by an executing of partialDFS.
                end %End while loop that checks whether to continue search
                    %from the vertex of B.
            end %End if statement checking whether vertex is free
        end %End of main DFS loop
       
    end %End main phase while loop
    
    %Return capacities in slightly different format.
    %Which corresponds to the format used in the Mapping code.
    capacity = zeros(2*n);
    capacity(n+1:2*n, 1:n) = capacityAB';
    capacity(1:n, n+1:2*n) = capacityBA';
end

% Utility function for Dijkstra's algorithm
function minIndex = DijkstraMinDistance(lv, shortestpathset, n)
    minValue = Inf;
    minIndex = -1; %placeholder value
    for v = 1:n*2
        if lv(v) < minValue && ~shortestpathset(v)
            minValue = lv(v);
            minIndex = v;
        end
    end
    
end

% Helper function that executes a single partial DFS from a free
% vertex of B. 
% vertex is the free vertex of B to initiate the search from.
% augmentingpathvertices is a column vector of length 2*n that describes
%   the current augmenting path's vertices, in order.
% AugPathVerIndex is the index of the last vertex of the augmenting path.
%   And is thus the number of vertices of the path and one larger than the
%   number of edges.
% vertexVisited is a 2 * n column vector, one slot per vertex v, where each
%   slot tracks the largest index of the explored neighbors of v.
%   vertexVisited must be retained over multiple executions of partialDFS
%   during a phase.
function [augmentingPathVertices, vertexVisited, AugPathVerIndex] = DFSUtil(vertex, augmentingPathVertices, n, slackAB, slackBA, AFree, vertexVisited, capacityAB, capacityBA)
    %Intialize the augmenting path as starting with a single vertex.
    AugPathVerIndex = 1;
    augmentingPathVertices(AugPathVerIndex) = vertex;
    
    %Keep going until we backtrack from the free vertex, meaning no
    %augmenting path was found.
    %If an augmenting path is found, we will break.
    while AugPathVerIndex > 0
        %Current vertex to explore from is the last vertex of the stack
        vertex = augmentingPathVertices(AugPathVerIndex);
        
        %If we reached a free vertex of A, return the path
        if vertex>n && AFree(vertex-n)
            return;
        end
        
        %Whether or not we need to backtrack due to not finding any
        %neighbors
        backtrack = true;
        
        %The first vertex in the adjacency of "vertex" to check next.
        range_var1 = vertexVisited(vertex)+1;
        %The last vertex in the adjacency of "vertex" to check.
        if vertex <= n
            %This is a vertex of B
            %We are exploring vertices of A as neighbors
            range_var2 = n*2;
        else
            %This is a vertex of A
            %We are exploring vertices of B as neighbors.
            range_var2 = n;
        end

        %Check all possible neighbors i of the current end of stack.
        for i = range_var1:range_var2
            %Update the farthest vertex explored from here.
            vertexVisited(vertex) = i;
            
            %Check if this edge is admissible. If so add it to the search
            %path
            if vertex <= n
                %vertex is type B
                a = i - n;
                if (slackBA(a,vertex) == 0) && (capacityBA(a,vertex) > 0)        
                    backtrack = false;
                    %Add vertex to path.
                    AugPathVerIndex = AugPathVerIndex + 1;
                    augmentingPathVertices(AugPathVerIndex) = i;  
                    break;
                end
            else
                %vertex is type A
                a = vertex - n;
                if (slackAB(i,a) == 0) && (capacityAB(i,a) > 0)        
                    backtrack = false;
                    %Add vertex to path
                    AugPathVerIndex = AugPathVerIndex + 1;
                    augmentingPathVertices(AugPathVerIndex) = i;  
                    break;
                end
            end
        end %End loop that checks all neighbors of vertex.
        if backtrack 
            %If no new edges outgoing from this vertex were found,
            %remove last vertex from path
            
            %This next line can theoretically be removed
            %but it makes things easier to debug.
            augmentingPathVertices(AugPathVerIndex) = 0;
            AugPathVerIndex = AugPathVerIndex - 1;
        end
    end
end












