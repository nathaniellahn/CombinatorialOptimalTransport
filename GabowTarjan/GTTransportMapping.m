% This function implements the mapping for the 1-feasibility
% implementation of the Gabow-Tarjan transportation algorithm, which
% acts on integral supplies, demands, and costs.
% This function calls the GTTransport function with the appropriate
% arguments. It also calls the GreedyMatch function to match the remaining
% supplies and demands not computed by GTTransport.


function [maxC, GTTransport_time, total_cost_transport, iterationCountTransport, APLengths, capacity_fulfilled] = GTTransportMapping(n, supplies, demands, CostInput, delta)
    
    % Mapping
    maxC = max(max(CostInput));
    CostNips=floor(4*CostInput/delta);
    newSupplies = floor((4*supplies*n*maxC)/delta);
    newDemands = ceil((4*demands*n*maxC)/delta);
    newDemands = newDemands';
    
    % Generating required cost matrix
    C = randi(n*2,n*2);
    for i = 1 : n
        for j = 1: n
            C(i,j) = 0;
        end
    end
    for i = n+1 : n*2
        for j = n+1 : n*2
            C(i,j) = 0;
        end
    end
    
    for i=1:n
        for j=n+1:n*2
            C(i,j) = CostNips(i,j-n);
        end
    end
    
    
    for i = n+1 : n*2
        for j = 1: n
            C(i,j) = C(j,i);
        end
    end
    
    
    timerValIn4 = tic;
    [iterationCountTransport, APLengths, capacity] = GTTransport(C, newSupplies, newDemands, n);
    
   
    actualSupplies = zeros(1,n);
    actualDemands = zeros(1,n);
    
    % Calculating matched supplies
    for i = 1:n
        actualDemands(i) =  sum(capacity(i+n,:));
        actualSupplies(i) = sum(capacity(:,i));
    end
    
    % Converting them to original inputs for comparison
    actualSupplies = (actualSupplies * delta)/(4*n*maxC);
    actualDemands = (actualDemands * delta)/(4*n*maxC);
    capacity(capacity~=0) = (capacity(capacity~=0) * (delta/(4*n*maxC)));
    
    
    % Identifying remaining supplies and demands
    remSupplies = supplies' - actualSupplies;
    remDemands = actualDemands - demands;

    % Pushing back extra demands
    for i=n+1:n*2
        j = 1;
        while remDemands(i-n)>0 && j<n+1 
           if capacity(i,j) >= remDemands(i-n)
               capacity(i,j) = capacity(i,j) - remDemands(i-n);
               remSupplies(j) = remSupplies(j) + remDemands(i-n);
               remDemands(i-n) = 0;
           elseif capacity(i,j) > 0
               remSupplies(j) = remSupplies(j) + capacity(i,j);
               remDemands(i-n) = remDemands(i-n) - capacity(i,j);
               capacity(i,j) = 0;
           end
           j = j+1;
        end
    end
    
    % Ensuring correct remaining demands after push back
    remDemands = abs(remDemands);
    
    % Computing greedy match for the remaining supplies and demands
    [greedyCapacity] = greedyMatch(n, remSupplies, remDemands, CostInput);
    timerValOut4 = toc(timerValIn4);
    GTTransport_time = timerValOut4;
    
    % Final edge capacities
    for i=1:n
        for j=n+1:n*2
            capacity(j,i) = capacity(j,i) + greedyCapacity(i,j-n);
        end
    end
   
   
   
    
    % Display solution
    %disp('Perfect Matching is');
    total_cost_transport = 0;
    capacity_fulfilled = 0;
    for i = 1:n
        for j = n+1:n*2
            if capacity(j,i) > 0
                capacity_fulfilled = capacity_fulfilled + capacity(j,i);
                total_cost_transport = total_cost_transport + (CostInput(i,j-n)*capacity(j,i));
                %fprintf('Vertex %d matched with %d with %d units and %d cost.\n',i,j, capacity(j,i), CostInput(i,j-n));
            end
        end
    end
    %fprintf('Total Cost of Matching is %f.\n',total_cost_transport);
    %fprintf('The number of iterations it took are: %d.\n',iterationCountTransport);
    %fprintf('The capacity fulfilled is: %f.\n',capacity_fulfilled);
    tolerance = 0.0000001;
    
    finalSuppliesUsed = zeros(n,1);
    finalDemandsUsed = zeros(n,1);
    for i = 1:n
        finalDemandsUsed(i) =  sum(capacity(i+n,:));
        finalSuppliesUsed(i) = sum(capacity(:,i));
    end
    
    
    residualS = abs(supplies - finalSuppliesUsed);
    residualD = abs(demands - finalDemandsUsed');
   
    %Verify solution is a valid transport plan
    if all(residualS <=  tolerance) ~= 1
       disp('Error: GT did not return a valid transport; a supply constraint was violated.')
    end

 
    if all(residualD <=  tolerance) ~= 1
       disp('Error2: GT did not return a valid transport; a demand constraint was violated.')
    end
    
    %Check total flow pushed is 1.
    assert(abs(capacity_fulfilled - 1) <= tolerance);

end