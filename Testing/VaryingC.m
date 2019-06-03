
% This file executes the set of tests that compare running times, 
% number of iterations, and length of augmenting paths of our algorithm
% for different values of C.

clc;

results = table;
n = 784;
global C;

%initialize parameters
delta = 0.1;
maxCosts = [100,200,300,400,500,600,700,800,900,1000];
runs = 10; %number of runs of the experiment

for k = 1:1:runs
     %% Verify solutions using LINPROG
    
    %random images
    disp('generate images');
    
    a = rand(n,1);
    b = rand(n,1);
    a = a/sum(a);
    b = b/sum(b);
    
    for maxCost = maxCosts
        disp(maxCost)
        %Generate a new cost matrix
        C = rand(n,n)*maxCost;
        
        %Scale so that the max cost becomes exactly the desired value
        C = C / max(max(C)) * maxCost;
        
        %% Run GTTransport
       
        [~,GTTransport_time,total_cost_transport,iterationCountTransport, APLengths,capacity_fulfilled] = GTTransportMapping(n, a, b', C, delta);

        %% Calculating Error
        errGT = total_cost_transport - lp_val;

        %% Saving results

        newRow = {n,k,maxCost,delta,total_cost_transport, GTTransport_time, iterationCountTransport, APLengths};
        results=[results;newRow];

    end
end

results.Properties.VariableNames = {'n', 'runNo', 'C', 'Delta','GTCost','GTTime','GTIter','GTAPLengths'};

%% Generate average results
avgResults = varfun(@mean,results,'InputVariables',{'n', 'Delta','GTCost','GTTime','GTIter','GTAPLengths'},'GroupingVariables',{'C'});


