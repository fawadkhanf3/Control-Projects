function [Ct] = win_always2(pwdyn, C0, rho, show_plot, verbose)
% win_always: compute set C ⊂ C0 such that
%   
%  C ⊂ Pre(C) + Ball(rho)
% 
% which is a sufficient condition for controlled invariance
%
% Reference: Rungger, M., & Tabuada, P. (2017). Computing Robust 
% Controlled Invariant Sets of Linear Systems. IEEE Trans. on 
% Automatic Control, dx.doi.org/10.1109/TAC.2017.2672859

  if nargin < 3
    rho = 0;
  end

  if nargin < 4
    show_plot = 0;
  end

  if nargin < 5
    verbose = 0;
  end

  
  if length(rho) == 1
    rho = rho * ones(pwdyn.n,1);
  end

  rho_ball = Polyhedron('A', [eye(pwdyn.n); -eye(pwdyn.n)], ...
                        'b', repmat(rho,2,1));

  C = Polyhedron('A', zeros(1,pwdyn.n), 'b', 1);
  Ct = C0;
  iter = 0;

  tic;

  if show_plot
    figure; clf
  end
  maxiter = 100;

%   while not (C-rho_ball <= Ct)
    while iter<maxiter
    C = Ct;

    if show_plot
      plot(C, 'alpha', 0.4); 
      drawnow
    end
    
    Cpre = pre(pwdyn,C, rho);
    Cpre = Cpre.merge;
    
    if isEmptyPolyUnion(Cpre)
      if verbose
        disp('returned empty')
      end
      Ct = Cpre;
      break;
    end
    
    Ct = intersect1(Cpre, C0);
%     Ct = minHRep(Ct);

%     cc = Ct.chebyCenter;
%     time=toc;
%     tic
        if isEmptySet(mldivide1(C, Ct))
			break;
		end
% toc
    iter = iter+1;
%     if verbose
%       disp(sprintf('iteration %d, %d ineqs, ball %f, time %f', ...
%                  iter, size(Ct.A,1), cc.r, time));
%     end
            disp(['Iteration #' num2str(iter) ' Complete.'])
            
  end

%   if verbose && ~isEmptySet(Ct)
%     disp(sprintf('finished with nonempty set after %d iterations!', iter))
%   end
