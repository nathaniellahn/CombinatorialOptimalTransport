% Note: This test file has been partially created using the testing code 
% from 
% https://github.com/chervud/AGD-vs-Sinkhorn
% and adapted for our purposes.

%This set of tests executes a comparison with the Sinkhorn algorithm.

%Import the compiled Java code for our algorithm.
%For the javaaddpath command, use the directory of the Java .class files.
clear java;
clc;
javaaddpath('..\GabowTarjanJavaCode\GTTransport\bin\');
import optimaltransport.*;

%Add all files of project to path
addpath(genpath('../../'));

load('mnist.mat', 'testX')

disp("Running in single threaded mode.")
maxNumCompThreads(1);

% example with set of written digits. Use digits.mat as input data
data = importdata('mnist.mat');
digits = testX';
results = table;

%size of an image
scale = 1; %multiplier for the number of pixels
m = 28*scale; %m x m - 2-dim grid, size of the support = m^2
%size of the support of the measures
n = m*m;

%Computes euclidean  distances.
C = computeDistanceMatrixGrid(m);

%Square the distances to get squared Euclidean distances.
C = C.*C;

%normalize cost matrix to set natural scale
%Ensures that C = 1;
C = C / max(max(C)); 
maxC = max(max(C));

%initialize parameters

deltas = [0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2];
runs = 10; %number of runs of the experiment

for k = 1:1:runs
    disp("Starting run #" + num2str(k));
    %load images from MNIST
    %generate image number
    
    i = ceil(rand*size(digits,2));
    j = ceil(rand*size(digits,2));
    
    imageIndex1 = i;
    imageIndex2 = j;
    aa = im2double(digits(:,i));
    bb = im2double(digits(:,j));
    aa = aa/sum(aa);
    bb = bb/sum(bb);
    aa = reshape(aa, m/scale, m/scale);
    aa = my_im_resize(scale,scale,aa);
    a = reshape(aa, m*m, 1);
    bb = reshape(bb, m/scale, m/scale);
    bb = my_im_resize(scale,scale,bb);
    b = reshape(bb, m*m, 1);
    b = b/sum(b);
    
    I = (a==0);
    a(I) =  0.000001;
    I = (b==0);
    b(I) =  0.000001;
    
    a = a/sum(a);
    b = b/sum(b);  
    
    %% Verify solutions using LINPROG
    
    disp('Computing Optimal Solution using LINPROG');
    lp_sol = -1;
    lp_val = -1;
    
    lp_time = tic;
    
    %If the optimal cost is not needed, comment this out as it takes
    %a long time to run.
    %[lp_sol,lp_val] = computeot_lp( C,a,b',n );
    lp_time = toc(lp_time);
 
    for delta = deltas
        %% Calculating their solution

        errSinkhorn = -1;
        sink2iters = -1;
        ot_val_sink2 = -1;
        
        %% Runs NIPS '17 version of Sinkhorn
        disp('start Sinkhorn');
        timerSink2 = tic;
        delta_prime = delta;
       
        eta = 4*log(n) / delta_prime;
        eps_prime = delta_prime / (8*maxC);
        A=exp(-1 * eta * C);
        A=A/sum(sum(A));
        [~, ot_val_sink2, sink2iters] = sinkhorn2(A, a, b', eps_prime, C);
        sink2iters = sink2iters - 1;
        timerSink2 = toc(timerSink2);
        errSink2 = ot_val_sink2 - lp_val;
        
        
        %% Run Greenkhorn
        ot_val_greenk = -1; 
        greenk_iters = -1;
        timerGreenk = -1;
        errGreenk = -1;
        
        disp("Start Greenkhorn")
        timerGreenk = tic();
        [~, ot_val_greenk, greenk_iters] = greenkhorn(A, a, b', eps_prime, C);
        greenk_iters = greenk_iters / n; %One "iter" for greenkhorn is n row column updates.
        timerGreenk = toc(timerGreenk);
        errGreenk = ot_val_greenk - lp_val;
        
        
        %% Run ICML '18 APDAGD
        AGDIters = -1;
        AGDTime = -1;
        
        disp("Start AGD")
        [AGDIters, AGDTime] = APDAGD(a, b, delta, C);
        
        
        %% Run GTTransport
        disp('Start GT Algorithm')
        
        gtSolver = optimaltransport.Mapping(n, a, b', C, delta);
        GTTransport_time = gtSolver.getTimeTaken();
        GTTransportMainRoutineTime = gtSolver.getMainRoutineTimeInSeconds();
       
        flow = gtSolver.getFlow();
        total_cost_transport = gtSolver.getTotalCost();
        iterationCountTransport = gtSolver.getIterations();
        APLengths = gtSolver.getAPLengths();

        %% Check to ensure that the solution is a valid transport plan.
        tolerance = 0.000000001;
        residualSupply = b';
        residualDemand = a;
        for i = 1:n
            for j = 1:n
                assert(flow(i,j) >= -tolerance);
                residualSupply(i) = residualSupply(i) - flow(i,j);
                residualDemand(i) = residualDemand(i) - flow(j,i);
            end
            assert(abs(residualSupply(i)) <= tolerance);
            assert(abs(residualDemand(i)) <= tolerance);
        end
 
        %% Calculating Error
        errGT = total_cost_transport - lp_val;
  
        %Verify that solution produced is sufficiently close to optimal
        %if the linprog code was run.
        if lp_val > 0
            assert(errGT <= delta);
        end
     
        %% Saving results
        newRow = {k, maxC, lp_val, lp_time, delta,timerSink2,errSink2,sink2iters,errGreenk,greenk_iters,timerGreenk,AGDIters, AGDTime,  GTTransport_time, iterationCountTransport, APLengths, errGT, imageIndex1, imageIndex2};
        results=[results;newRow];
        
    end
end


results.Properties.VariableNames = {'run', 'C','LINPROGCost','LINPROGTime','delta','SinkTime', 'SinkError','sinkiters','errGreenk','greenk_iters','timerGreenk', 'AGDIters', 'AGDTime', 'GTTime', 'GTIter','GTAPLengths' 'ErrorGT', 'ImageOneIndex', 'ImageTwoIndex'};

%% Generate average results
avgResults = varfun(@mean,results,'InputVariables',{'C','LINPROGCost','LINPROGTime', 'SinkTime', 'SinkError','sinkiters', 'errGreenk','greenk_iters','timerGreenk','AGDIters', 'AGDTime','GTTime', 'GTIter','GTAPLengths' 'ErrorGT'},'GroupingVariables',{'delta'});
disp(avgResults);



