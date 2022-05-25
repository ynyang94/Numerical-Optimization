function [x,y,z,B,flag,iter] = mysimplex(A,b,c,x,B,maxiter,varargin)
    % function [x,y,z,B,flag,iter] = mysimplex(A,b,c,x,B,maxiter)
    % function [x,y,z,B,flag,iter] = mysimplex(A,b,c,x,B,maxiter,debug)
    % 
    % MYSIMPLEX Solves the LP
    %
    %       min  c'*x 
    %       s.t. Ax=b, x>=0.
    %
    % using the simplex method
    %
    % inputs:
    %  A	mxn matrix (m<n)
    %  b	m vector (RHS of constraints)
    %  c	n vector (objective)
    %  x	n vector basic feasible point
    %  B	m vector of indices corresponding to basis of x
    %  maxiter  maximum number of iterations (maxiter = inf is allowed)
    %  debug [optional]  print intermediate results if debug == 1
    %		   
    % outputs:
    %  x	n vector optimal solution to the LP   (only if flag==0)
    %  y	m vector optimal Lagrange multipliers (only if flag==0)
    %  z	n vector optimal slacks               (only if flag==0)
    %  B	m vector of indices corresponding to basis of x
    %  flag	0 if optimal solution was found
    %       1 if problem is unbounded
    %       2 if the number of iterations is exceeded
    %  iter number of simplex iterations executed
    %
    
    
 
     if (nargin == 6)
         debug=0;
     else
         debug = varargin{1};
     end
         
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % sanity checks
      [m,n]=size(A); 
     if (m>=n) error('must have m<n'); end
     if (m~=size(b,1) | 1~=size(b,2)) error('b must be an mx1 vector'); end
     if (n~=size(c,1) | 1~=size(c,2)) error('c must be an nx1 vector'); end
     if (rank(A)<m) error('A must be full rank'); end
    
     if (length(B)~=m)==1
         error('start bfp does not have exactly m non zero entries');
     end
     % some checks to see if (x0,B) is really a bfp
     % check feasability of x0
     if (any(x<-10*eps)) 
         error('x not feasible: it has negative entries!!'); 
     end
     if (norm(A*x-b,inf)>10*eps) 
         error('x not feasible: Ax0 != b '); 
     end
     % check all components outside basis are zero
     if ( sum(x(setdiff(1:n,B)))>10*eps )
         error('not all components of x outside basis set are zero');
     end
     % check A(:,B) is invertible
     if (rank(A(:,B))<m)
         error('the columns of A corresp. to this bfp are not lin. indep');
     end
     
     %implement simplex method
     iter=0;
     z=zeros(n,1);
     %set initial N={1,2,..6}\B
     N=setdiff(1:6,B);
     %pirnt out initial value.
     fprintf("Iterations:%d\n",iter)
     %compute objective value.
     obj=c'*x;
     fprintf("Objective function:%.2f\n",obj)
     disp(sprintf("Basis: %d, %d, %d\n", B))
     disp(sprintf("Vertex: %.2f, %.2f, %.2f,%.2f,%.2f,%.2f\n", x))
     
     %implement simplex method in Notes4-10, Algorithm 4.43
     while iter<maxiter
         
         %divide matrix by index
         c_b=c(B);
         c_n=c(N);
         A_b=A(:,B);
         A_n=A(:,N);
         x_b=x(B);
         x_n=x(N);
         
         %compute y and z_n
         y=inv(A_b)'*c_b;
         disp(sprintf('Lagrange Multiplier y:%.2f,%.2f,%.2f\n', y))
         
         %construct z_n
         z_n=c_n-transpose(A_n)*y;
         %construct z.
         z(B)=0;
         z(N)=z_n;    
         %print out z.
         disp(sprintf('Lagrange Multiplier z:%.2f, %.2f, %.2f,%.2f,%.2f,%.2f\n', z))
         %check optimality condition
         %set the initial flag be 0.
         flag=0;
         if all(z_n >=0)
             %optimality condition achieves
 
             flag=0;
             %stop iteration
             break
         
         else
             %find index j0 with z_n<0;
             [zj0,j0]=min(z_n);
             
             %compute d_b
             d=inv(A_b)*A(:,N(j0));
             D(B)=d;
             D(N)=0;
             D(j0)=-1;
             disp(sprintf('direction d:%.2f, %.2f, %.2f,%.2f,%.2f,%.2f\n', D))
             i_index=find(d>0);
             d_b=d(find(d>0));
             if all(d_b<=0)==1
                 flag=1;
                 %stop iteration
                 break
                 %this case, the optimal value is unbounded.
                 
             else
                 %find t.
                 update=x_b(i_index)./d_b;
                 %compute the largest value t can achieve
                 t=min(update);
             end
         end
         %set new solution.
         x(B)=x(B)-(t*d);
         %!
         x(N(j0))=t;
         %set new index set
         B=find(x>0);
         N=find(x==0);
         %construct z
         %update iteration.
         iter=iter+1;
         obj=c'*x;
         %print out some solution
         fprintf("Iteration:%d\n",iter)
         disp(sprintf('solution x:%.2f, %.2f, %.2f,%.2f,%.2f,%.2f\n', x))
         fprintf('objective function is:%d\n',obj);
         disp(sprintf("Basis: %d, %d, %d\n", B))
         
     end
     if (iter==maxiter)
         flag=2;
     end
end
