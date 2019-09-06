function [iter,time,A,B]=Sinkhorn(r,c,epsilon)

%An implementation of Sinkhorn's algorithm 
%from [Altshuler et al, 2017], Alg.3
%Code was retrieved from https://github.com/chervud/AGD-vs-Sinkhorn
%with some minor changes:

%1.) The 20,000 iteration cap was removed. For the tests executing in AGD,
%it seemed like 20,000 iterations were sufficient, but we did not want to
%make this assumption in general. Some tests with very low error may
%require 20,000 iterations to ensure that the proper error value is
%reached.
%2.) We added a final rounding step. Before, the AGD implementation did not
%include this final rounding step. This final rounding step is necessary 
%to know the computed tranpsort cost and error, but has very little effect
%on the running time.


%B is the transport plan prior to the final rounding step.
%A is the final transport plan after the rounding step (the one that
%should be used in nearly all cases.).


% r -- vector of first measure
% c -- vector of second measure

tic;


n = size(r,1);

global C; %cost matrix

gamma = epsilon/4/log(n); %regularization parameter

A = exp(-C/gamma);
A_0 = A;
A = A/sum(sum(A));

max_el = max(max(C));

k=0;
x = zeros(n,1);
y = zeros(n,1);

while (norm(r - sum(A,2),1) + norm(c - sum(A,1)',1) > epsilon/8/max_el)
    k = k+1;
    if mod(k,2) == 1
        x = x + log(r./sum(A,2));
    else
        y = y + log(c./(sum(A,1)'));
    end
    C_new = -C/gamma+x*ones(1,n)+ones(n,1)*y';
    A = exp(C_new);
end


%Round off the supplies and demands to satisfy the feasibility constraints
%Due to the while loop stopping condition, the cost error incurred is
%guaranteed to be acceptable.
%Note, this was added to the Sinkhorn code from the AGD implementation,
%because it was mainly interested in the number of iterations,
%not the final solution; they did not execute the final rounding step for
%Sinkhorn. 

%Store matrix prior to rounding: this is what the AGD implementation
%returned
B = A;

%Round the matrix to polytope to produce a valid transport.
A = round_matrix(A,r,c);
time = toc;

%str = ['average time per iteration ',num2str(time/k)];
%disp(str); %print current iteration number
iter = k;

end