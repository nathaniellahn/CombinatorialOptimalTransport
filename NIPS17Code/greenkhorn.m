%% 
% NOTE: This implementation was retrieved from 
%   https://github.com/JasonAltschuler/OptimalTransportNIPS17
% We made a few modifications, all of which were meant to address the 
% stopping condition. In particular, the original implementation 
% ran a fixed number of iterations, which does not guarantee an error
% of delta or less. Instead, we added the stopping condition given
% in their paper.
%
% The following comments are from their code base, modified slightly
% to describe new parameter changes.
% Implementation of our Greenkhorn algorithm for matrix scaling.
% Specifically, each iteration, choose the row or column that gives
% maximum "gain" (w.r.t. our gain function), and then normalize P to have
% total mass 1. See paper for further details.
%
% Input parameters:
%  -- A:  positive square matrix to project onto U_{r,c} transport polytope (dims: nxn)
%  -- r:  desired row sums (marginals)         (dims: nx1)
%  -- c:  desired column sums (marginals)      (dims: 1xn)
%  -- eps_prime: the minimum distance to polytope necessary to 
%       ensure a good error.
%  -- C:  cost matrix for OT
% 
% Output:
%  -- P:   final scaled matrix
%  -- otval:  optimal transport value
%  -- t: number of iterations.

function [P, otval, t] = greenkhorn(A,r,c,eps_prime,C)
    
    addpath(genpath('algorithms/'));
    P = A;

    % compute full row and column marginals once
    r_P = sum(P,2);
    c_P = sum(P,1);

    % compute gains for each row and column
    r_gain = r_P - r + r.*log(r./r_P);
    c_gain = c_P - c + c.*log(c./c_P);

    
    err = norm(r_P-r,1)+norm(c_P-c,1);

    t = 0;
    while err > eps_prime
        
        t = t + 1;
        % find row or column with maximum gain
        [r_gain_max, i] = max(r_gain);
        [c_gain_max, j] = max(c_gain);

        if r_gain_max > c_gain_max        
            % update row i
            scaling = r(i)/r_P(i);
            old_row = P(i,:);
            new_row = old_row*scaling;
            P(i,:)  = new_row;

            % renormalize (can also be done implicitly if one wants to optimize)
            P = P/sum(sum(P));

            % compute full row and column marginals
            r_P = sum(P,2);
            c_P = sum(P,1);

            % compute gains for each row and column
            r_gain = r_P - r + r.*log(r./r_P);
            c_gain = c_P - c + c.*log(c./c_P);

            % % tricks to speed up computation if we are not renormalizing
            % % matrix each time
            %         % update row and column sums in O(n) time
            %         r_P(i)  = r(i);
            %         c_P     = c_P - old_row + new_row;
            %
            %         % update row and column gains in O(n) time
            %         r_gain(i) = 0;
            %         c_gain    = c_P - c + c.*log(c./c_P);

            err = norm(r_P-r,1)+norm(c_P-c,1);        
        else
            % update column j
            scaling = c(j)/c_P(j);
            old_col = P(:,j);
            new_col = old_col*scaling;
            P(:,j)  = new_col;

            % renormalize (can also be done implicitly if one wants to optimize)
            P = P/sum(sum(P));

            % compute full row and column marginals
            r_P = sum(P,2);
            c_P = sum(P,1);

            % compute gains for each row and column
            r_gain = r_P - r + r.*log(r./r_P);
            c_gain = c_P - c + c.*log(c./c_P);

            % % tricks to speed up computation if we are not renormalizing
            % % matrix each time
            %       % update row and column sums in O(n) time
            %       c_P(j)  = c(j);
            %       r_P     = r_P - old_col + new_col;
            %
            %       % update row and column gains in O(n) time
            %       c_gain(j) = 0;
            %       r_gain = r_P - r + r.*log(r./r_P);

            err = norm(r_P-r,1)+norm(c_P-c,1);
        end
        if mod(t, 1000) == 0
            %disp("Error: " + num2str(err) + ", Goal:" + num2str(eps_prime));
        end
    end
    otval =  frobinnerproduct(round_transpoly(P,r,c),C);
end