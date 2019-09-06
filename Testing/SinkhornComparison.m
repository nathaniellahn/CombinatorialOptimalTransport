% Note: This test file has been partially created using the testing code 
% from 
% https://github.com/chervud/AGD-vs-Sinkhorn
% and adapted for our purposes.

%This set of tests executes a comparison with the Sinkhorn algorithm.

clc;

% example with set of written digits. Use digits.mat as input data
data = importdata('mnist.mat');
digits = testX';
results = table;
%size of an image
scale = 1; %multiplier for the number of pixels
m = 28*scale; %m x m - 2-dim grid, size of the support = m^2
%size of the support of the measures
n = m*m;

% calculate transport cost matrix
global C;
disp('computeDistanceMatrixGrid');

%Computes eudclidean  distances.
C = computeDistanceMatrixGrid(m);

%Square the distances to get squared Euclidean distances.
C = C.*C;

%normalize cost matrix to set natural scale
C = C / median(C(:)); 

%initialize parameters
disp(max(max(C)))

%deltas = [0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2];
%deltas = [0.2,0.4,0.6,0.8,1.0]
deltas = [0.05, 0.1, 0.15,0.2, 0.25];

runs = 10; %number of runs of the experiment

for k = 1:1:runs
    %load images from MNIST
    %generate image number
    i = ceil(rand*size(digits,2));
    j = ceil(rand*size(digits,2));
 
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
    
    disp('LINPROG SOLUTION');
    timerValIn3 = tic;
   
    lp_sol = -1;
    lp_val = -1;
    
    %If the optimal cost is not needed, comment this out as it takes
    %a long time to run.
    [lp_sol,lp_val] = computeot_lp( C,a,b',n );
   
    timerValOut3 = toc(timerValIn3);
    lp_time = timerValOut3;
   
    
    
    for delta = deltas
        errGTFlag = 0;
        %%
        %% Run Sinkhorn
        disp('start Sinkhorn');
        
        %[iter,time,A,~] = Sinkhorn(a,b,delta);
        
        %Run this code if we want to give sinkhorn a more lenient delta.
        %Such that the actual errors generated match ours.
        [iter,time,A,~] = Sinkhorn(a,b,3.5*delta);
       
        %% Calculating their solution
        
        solSinkhorn = 0;
        for i=1:n
            for j=1:n
                solSinkhorn = solSinkhorn+(A(i,j)*C(i,j));
            end
        end
        errSinkhorn = solSinkhorn - lp_val;
        %check to ensure that the solution is a valid transport plan.
        tolerance = 0.0000001;
        residualR = abs(a - sum(A,2));

        if all(residualR <=  tolerance) ~= 1
           disp('Error:R Sinkhorn did not return a valid transport')
        end
        
        residualC = abs(b - sum(A,1)');
       
        if all(residualC <=  tolerance) ~= 1
           disp('Error:C Sinkhorn did not return a valid transport')
        end
    
        
        %% Run GTTransport
        disp('GT Transport SOLUTION');
        
        
        %The code to run if your only goal is to get a transport with error
        %at most \delta.
        [~,GTTransport_time,total_cost_transport,iterationCountTransport, APLengths,capacity_fulfilled] = GTTransportMapping(n, a, b', C, delta);
        
        %The code to execute for if we want to match roughly match
        %Sinkhorn's error
        %same as above except uses delta/6 as input.
        %[maxCost,GTTransport_time,total_cost_transport,iterationCountTransport, APLengths,capacity_fulfilled] = GTTransportMapping(n, a, b', C, delta/6);
        
        
        
        %% Calculating Error
        errGT = total_cost_transport - lp_val;
        
        if errGT > delta
            errGTFlag = 1;
        end
        
        %% Saving results
        maxCost = max(max(C));
        newRow = {k, maxCost, lp_val, lp_time, delta,time, iter, errSinkhorn, GTTransport_time, iterationCountTransport, APLengths, errGT, i, j};
        results=[results;newRow];
        
    end
end


results.Properties.VariableNames = {'run', 'C','LINPROGCost','LINPROGTime','delta','SinkhornTime','SinkhornIter','ErrorSinkhorn', 'GTTime', 'GTIter','GTAPLengths' 'ErrorGT','ImageOneIndex', 'ImageTwoIndex'};

%% Generate average results
avgResults = varfun(@mean,results,'InputVariables',{'C','LINPROGCost','LINPROGTime','SinkhornTime','SinkhornIter','ErrorSinkhorn', 'GTTime', 'GTIter','GTAPLengths' 'ErrorGT'},'GroupingVariables',{'delta'});


