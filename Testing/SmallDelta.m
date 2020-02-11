% Note: This test file has been partially created using the testing code 
% from 
% https://github.com/chervud/AGD-vs-Sinkhorn
% and adapted for our purposes.

%This set of tests executes a comparison with the Sinkhorn algorithm.

%Import the compiled Java code for our algorithm.
%For the javaaddpath command, use the directory of the Java .class files.
clear java;
clc;

%Add all files of project to path
addpath(genpath('../'));

javaaddpath('..\GabowTarjanJavaCode\GTTransport\bin\');
import optimaltransport.*;

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

deltas = [0.1, 0.01, 0.001, 0.0001];

runs = 100; %number of runs of the experiment

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

    for delta = deltas
        %% Run GTTransport
        disp('Start GT Algorithm')
        
        gtSolver = optimaltransport.Mapping(n, a, b', C, delta);
        GTTransport_time = gtSolver.getTimeTaken();
        GTTransportMainRoutineTime = gtSolver.getMainRoutineTimeInSeconds();
       
        flow = gtSolver.getFlow();
        total_cost_transport = gtSolver.getTotalCost();
        iterationCountTransport = gtSolver.getIterations();
        APLengths = gtSolver.getAPLengths();
        augmentTime = gtSolver.getTimeTakenAugment();
        
        %augmentTime is 0 if the timing code is commented out in the Java 
        %implementation.
        %Precise timing could negatively affect performance.
        %If the time taken from augmentations is important,
        %then uncomment the timing calls in the Java code.
        if augmentTime == 0
            augmentTime = -1;
        end

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
        
        %% Saving results
        newRow = {k, maxC, delta, GTTransport_time, GTTransportMainRoutineTime, iterationCountTransport, APLengths, augmentTime, imageIndex1, imageIndex2};
        results=[results;newRow];
        
    end
end


results.Properties.VariableNames = {'run', 'C','delta','GTTime', 'GTTransportMainRoutineTime', 'GTIter','GTAPLengths', 'AugmentTime','ImageOneIndex', 'ImageTwoIndex'};

%% Generate average results
avgResults = varfun(@mean,results,'InputVariables',{'GTTime', 'GTTransportMainRoutineTime', 'GTIter','GTAPLengths', 'AugmentTime' },'GroupingVariables',{'delta'});
disp(avgResults);



