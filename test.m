%%test case
%set initial value of A,b,c
A = [-2, 1, 1, 1, 0, 0;
     -1, 2, 0, 0, 1, 0;
      1, 0, 0, 0, 0, 1;];
b = [2; 7; 3];
c = [-1; -2; 0; 0; 0; 0];
x = [0; 0; 0; 2; 7; 3];
B = [ 4; 5;6;];
%set maximal iteration.
maxiter = 10;
%solve LP problem by simplex method 
[x,y,z,B,flag,iter] = mysimplex(A,b,c,x,B,maxiter,1);
%print out solution.
fprintf('mysimplex terminated with \n');
fprintf('flag=%d\n',flag);
fprintf('iter=%d\n',iter);
fprintf('x=['); fprintf('%d ',x); fprintf(']\n');
fprintf('y=['); fprintf('%d ',y); fprintf(']\n');
fprintf('z=['); fprintf('%d ',z); fprintf(']\n');  


% %test the answer is 3,5,3
% A = [-2, 1, 1;
%     -1, 2, 0;
%     1, 0, 0;];
% c=[-1;-1;0];
% b=[2;7;3];
% x=linprog(c,A,b)




% function [x,y,z,B,flag,iter] = mysimplex(A,b,c,x,B,maxiter,varargin)
%     % function [x,y,z,B,flag,iter] = mysimplex(A,b,c,x,B,maxiter)
%     % function [x,y,z,B,flag,iter] = mysimplex(A,b,c,x,B,maxiter,debug)
%     % 
%     % MYSIMPLEX Solves the LP
%     %
%     %       min  c'*x 
%     %       s.t. Ax=b, x>=0.
%     %
%     % using the simplex method
%     %
%     % inputs:
%     %  A	mxn matrix (m<n)
%     %  b	m vector (RHS of constraints)
%     %  c	n vector (objective)
%     %  x	n vector basic feasible point
%     %  B	m vector of indices corresponding to basis of x
%     %  maxiter  maximum number of iterations (maxiter = inf is allowed)
%     %  debug [optional]  print intermediate results if debug == 1
%     %		   
%     % outputs:
%     %  x	n vector optimal solution to the LP   (only if flag==0)
%     %  y	m vector optimal Lagrange multipliers (only if flag==0)
%     %  z	n vector optimal slacks               (only if flag==0)
%     %  B	m vector of indices corresponding to basis of x
%     %  flag	0 if optimal solution was found
%     %       1 if problem is unbounded
%     %       2 if the number of iterations is exceeded
%     %  iter number of simplex iterations executed
%     %
%     
%     
%  
%      if (nargin == 6)
%          debug=0;
%      else
%          debug = varargin{1};
%      end
%          
%      
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      % sanity checks
%      [m,n]=size(A); 
%      if (m>=n) error('must have m<n'); end
%      if (m~=size(b,1) | 1~=size(b,2)) error('b must be an mx1 vector'); end
%      if (n~=size(c,1) | 1~=size(c,2)) error('c must be an nx1 vector'); end
%      if (rank(A)<m) error('A must be full rank'); end
%     
%      if (length(B)~=m)
%          error('start bfp does not have exactly m non zero entries');
%      end
%      % some checks to see if (x0,B) is really a bfp
%      % check feasability of x0
%      if (any(x<-10*eps)) 
%          error('x not feasible: it has negative entries!!'); 
%      end
%      if (norm(A*x-b,inf)>10*eps) 
%          error('x not feasible: Ax0 != b '); 
%      end
%      % check all components outside basis are zero
%      if ( sum(x(setdiff(1:n,B)))>10*eps )
%          error('not all components of x outside basis set are zero');
%      end
%      % check A(:,B) is invertible
%      if (rank(A(:,B))<m)
%          error('the columns of A corresp. to this bfp are not lin. indep');
%      end
%     
%  
% end
