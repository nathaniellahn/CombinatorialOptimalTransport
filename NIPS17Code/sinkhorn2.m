%% 
% Implementation of classical Sinkhorn algorithm for matrix scaling.
% Each iteration simply alternately updates (projects) all rows or
% all columns to have correct marginals.
%
% NOTE: This implementation was retrieved from 
%   https://github.com/JasonAltschuler/OptimalTransportNIPS17
% We made a few modifications, all of which were meant to address the 
% stopping condition. In particular, the original implementation 
% ran a fixed number of iterations, which does not guarantee an error
% of delta or less. Instead, we added the stopping condition given
% in their paper.
% 
% Input parameters:
%  -- A:  positive square matrix to project onto U_{r,c} transport polytope (dims: nxn)
%  -- r:  desired row sums (marginals)         (dims: nx1)
%  -- c:  desired column sums (marginals)      (dims: 1xn)
%  -- eps_prime: Minimum distance to polytope necessary to stop
%  -- compute_otvals: flag whether to compute otvals (slow but used in some plots)
%  -- C:  cost matrix for OT
%
% Output:
%  -- P:   final scaled matrix
%  -- err: sum of row and column violations at each iteration
%  -- ot:  values of optimal transport of matrix iterates

function [P, otval, t] = sinkhorn2(A,r,c,eps_prime,C)
P = A;
r_P = sum(P,2);
c_P = sum(P,1);
err = norm(r_P-r,1)+norm(c_P-c,1);

%We changed the stopping condition from the original implementation
%Because the original ran some fixed number of iterations, while
%our interest is producing an error of at most delta.
t = 1;
while err > eps_prime
    if mod(t,2)==1
        % rescale rows
        r_P = sum(P,2);
        P   = bsxfun(@times,P,r./r_P);
        
        r_P = sum(P,2);
        c_P = sum(P,1);
        err = norm(r_P-r,1)+norm(c_P-c,1);
    else
        % rescale columns
        c_P = sum(P,1);
        P   = bsxfun(@times,P,c./c_P);
        
        r_P = sum(P,2);
        c_P = sum(P,1);
        err = norm(r_P-r,1)+norm(c_P-c,1);
    end
    t = t + 1;
    %disp("Error: " + num2str(err) + ", Goal:" + num2str(eps_prime));
end
otval = frobinnerproduct(round_transpoly(P,r,c),C);
end