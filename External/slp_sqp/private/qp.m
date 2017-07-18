function [s,u,status,steplimitactive]=qp(H,fp,A,b,vlb,vub,x0,neq,nomsg,max_s)
% QP converts old optimization toolbox qp calling sequence to new quadprog
% calling sequence in version 2 and later of the optimization toolbox.
%
% *If using version 1 optimization toolbox, remove this function.*

%--Modifications
%   5/12/06 - Structure options
%    4/4/07 - Increase MaxIter
%    5/6/13 - If A is null, then neq=0
%   9/30/15 - Default quadprog alogorithm instead of deprecated active-set
%  10/15/15 - default quadprog fails for nonconvex--revert to   active-set
%  11/16/15 - default quadprog no "s" for infeasible--revert to active-set
%  11/25/16 - active-set algorithm removed for R2106b

if nargin < 10
    max_s = [];
end

[nc,nx]=size(A);
neq=min(neq,nc);
MaxIter=max(200,10*min(nx,100));
% Rows of A are partitioned into equality and inequality constraints: p1, p2
p1=1:neq;
p2=neq+1:nc;
% No displayed output (nomsg = -1)
msg={'off','final','iter'};
%%%%%MOSEK%%%%%%
%    path('c:\mosek\4\toolbox\r14sp3',path)
%   options=optimset('Display',msg{nomsg+2},'MaxIter',5000);
%%%%%NORMAL%%%%%%
options=optimset('quadprog');
if verLessThan('optim','7.0')
   options=optimset(options,'Display',msg{nomsg+2},'LargeScale','off',...
                    'Algorithm','active-set','MaxIter',MaxIter);
else
   options=optimset(options,'Display',msg{nomsg+2},'MaxIter',MaxIter);
end

% Set additional lower and upper bounds to reflect limits in allowed step size for
% each variable
nx = numel(fp);
if isempty(max_s)
    A_slim = zeros(0,nx);
    b_slim = zeros(0,1);
else
    [~,~,V] = svd(H);
    A_slim = [-V.'; V.'];
    b_slim = repmat(max_s, 2*nx, 1);
end

if isempty(H)
   [s,f,exitflag,output,LAMBDA]=linprog(fp,[A(p2,:);A_slim],[b(p2);b_slim],A(p1,:),b(p1),vlb,vub,x0,options);
else
   [s,f,exitflag,output,LAMBDA]=quadprog(H,fp,[A(p2,:);A_slim],[b(p2);b_slim],A(p1,:),b(p1),vlb,vub,x0,options);
   if exitflag<0 && verLessThan('optim','7.5')
      warning('off','optim:quadprog:WillBeRemoved');
      options=optimset(options,'Algorithm','active-set');
      [s,f,exitflag,output,LAMBDA]=quadprog(H,fp,[A(p2,:);A_slim],[b(p2);b_slim],...
         A(p1,:),b(p1),vlb,vub,x0,options);
   elseif exitflag==-6 || exitflag==0
      if isempty(s), s=zeros(size(fp)); end
      HessianFcn = @(x,lambda) H;
      options = optimoptions(@fmincon,'SpecifyObjectiveGradient', true,...
         'HessianFcn',HessianFcn,...
         'Display','off');
      [s,f,exitflag,output,LAMBDA]=fmincon(@qobj,s,[A(p2,:);A_slim],[b(p2);b_slim],...
         A(p1,:),b(p1),vlb,vub,[],options);
   end
end
if exitflag>0
   status='ok';
elseif exitflag==0
   status='maximum number of iterations exceeded';
elseif exitflag<0
   status='infeasible, unbounded, or unconverged';
end
u=[];
if exitflag>=0
   if ~isempty(LAMBDA.eqlin),   u=[u;LAMBDA.eqlin(:)]; end
   if ~isempty(LAMBDA.ineqlin), u_ineq=LAMBDA.ineqlin([p1 p2]); u=[u;u_ineq(:)]; end
   if ~isempty(vlb),            u=[u;LAMBDA.lower(1:length(vlb))]; end
   if ~isempty(vub),            u=[u;LAMBDA.upper(1:length(vub))]; end
end
steplimitactive = ~isempty(max_s) && any(abs(V.'*s) >= 0.999*max_s);
%%Added -1 conditions to handle MOSEK exitflags
%  if exitflag==-1
%     status='ok';
%  end
% %%%%%MOSEK%%%%%%%
%  rmpath('c:\mosek\4\toolbox\r14sp3')
%  u = [LAMBDA.eqlin(:); LAMBDA.ineqlin(:); ...
%       LAMBDA.lower(1:length(vlb)); LAMBDA.upper(1:length(vub))];

%% Quadratic objective function for QP search direction subproblem
   function [f,gradf]=qobj(x)
      f     = fp'*x + x'*H*x/2;
      gradf = fp    +    H*x;
   end
end