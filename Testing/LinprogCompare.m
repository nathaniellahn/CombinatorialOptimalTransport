% Note: This test file has been partially created using the testing code 
% from 
% https://github.com/chervud/AGD-vs-Sinkhorn
% and adapted for our purposes.

% This test executes the GT Algorithm on small delta values, and compares
% the resulting transport costs and running times to linprog.
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
C = computeDistanceMatrixGrid(m);
disp('C.*C');
C = C.*C;
disp('C / median(C(:))');
C = C / median(C(:)); %normalize cost matrix to set natural scale for \gamma

%initialize parameters

epsilon = [0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035, 0.004, 0.0045, 0.005];

runs = 10; %number of runs of the experiment

for k = 1:1:runs
    %% Select 2 random images from MNIST data set.
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
    
    %Record time taken by linprog
    timerValIn3 = tic;
    
    [lp_sol,lp_val] = computeot_lp( C,a,b',n );
    timerValOut3 = toc(timerValIn3);
    lp_time = timerValOut3;
   
    
    
    for eps = epsilon
        errGTFlag = 0;
        
        %% Run GTTransport
        disp('GT Transport SOLUTION');
       
        delta = eps;
        [~,GTTransport_time,total_cost_transport,iterationCountTransport, APLengths,capacity_fulfilled] = GTTransportMapping(n, a, b', C, delta);
        
        %% Calculating Error
        errGT = total_cost_transport - lp_val;
        
        if errGT > eps
            errGTFlag = 1;
        end
        
        %% Saving results
        
        maxCost = max(max(C))
        newRow = {k, maxCost, lp_val, lp_time, eps, GTTransport_time, iterationCountTransport, APLengths, errGT, i, j};
        results=[results;newRow];
        
    end
end

results.Properties.VariableNames = {'run', 'C','LINPROGCost','LINPROGTime','delta', 'GTTime', 'GTIter','GTAPLengths' 'ErrorGT', 'ImageOneIndex', 'ImageTwoIndex'};

%% Generate average results
avgResults = varfun(@mean,results,'InputVariables',{'C','LINPROGCost','LINPROGTime', 'GTTime', 'GTIter','GTAPLengths' 'ErrorGT'},'GroupingVariables',{'delta'});


