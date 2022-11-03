function [ X0 ] = pre_proj(dyn, X, rho)
  % Compute set
  % {x : ∀ p ∃ u ∀ d, x(t+1) + Ball(rho) ⊆ X}
  % using normal intersection/projection method
  %
  % Reference: Petter Nilsson Ph.D. thesis (2017), Theorem 3.4

  if ~isa(dyn, 'Dyn')
    error('dyn must be an instance of Dyn');
  end

  if nargin < 3
    rho = zeros(dyn.nx,1);%0;
  end

  if length(rho) == 1
    rho = rho * ones(dyn.nx,1);
  end

  % check if X is empty set
  if(isEmptySet(X))
      X0 = X;
      return;
  end

  Xb = X - Polyhedron('A', [eye(dyn.nx); -eye(dyn.nx)], 'b', repmat(rho,2,1));
  Xb.minHRep; % not sure it is necessary or not. 
  
  X0_A = zeros(0, dyn.nx);
  X0_b = ones(0,1);

  PV = dyn.P.V;
  DV = dyn.D.V;

  N_P = max(1, size(PV,1));       % number of p vertices
  N_V = max(1, length(dyn.XV_V)); % number of v vertices

  proj = cell(N_P * N_V, 1);

  % For each (p,v) combination
  for iter=1:N_P*N_V
    ip = 1+mod(iter-1, N_P);    % idx p
    iv = 1+floor((iter-1)/N_P); % idx v

    A_mat_p = dyn.A;
    F_mat_p = dyn.F;

    if dyn.nv > 0
      vext_x = dyn.XV_V{iv}(:, 1:dyn.nx);
      vext_v = dyn.XV_V{iv}(:, dyn.nx+dyn.np+1:end);
      A_mat_p = A_mat_p + dyn.Ev * vext_x;
      F_mat_p = F_mat_p + dyn.Ev * vext_v;
      if dyn.np > 0
        vext_p = dyn.XV_V{iv}(:, dyn.nx+1:dyn.nx+dyn.np);
        F_mat_p = F_mat_p + dyn.Ev * vext_p*PV(ip, :)';
      end
    end

    for jp=1:dyn.np
      A_mat_p = A_mat_p + dyn.Ap{jp} * PV(ip, jp);
      F_mat_p = F_mat_p + dyn.Fp{jp} * PV(ip, jp);
    end

    Xd_A = zeros(0, dyn.nx+dyn.nu);
    Xd_b = zeros(0, 1);

    % For each d
    for id=1:max(1, size(DV,1))
      A_mat_pd = A_mat_p;
      F_mat_pd = F_mat_p;
      for jd=1:dyn.nd
        A_mat_pd = A_mat_pd + dyn.Ad{jd} * DV(id, jd);
        F_mat_pd = F_mat_pd + dyn.Fd{jd} * DV(id, jd);
      end
      if dyn.nw > 0
        % For each w
        for iw=1:length(dyn.XW_V)
          wext_x = dyn.XW_V{iw}(:, 1:dyn.nx);              % Hdx
          wext_p = dyn.XW_V{iv}(:, dyn.nx+1:dyn.nx+dyn.np);
%           wext_w = dyn.XW_V{iw}(:, dyn.nx+dyn.np+1:end);
          wext_w = dyn.XW_V{iw}(:, dyn.nx+1);              % hd
          wext_u = dyn.XW_V{iw}(:, dyn.nx+dyn.np+1+1:end); % Hdu
%           Xd_A = [Xd_A; 
%                   Xb.A*[A_mat_pd+dyn.Ew*wext_x dyn.B]];
              
          Xd_A = [Xd_A; 
                  Xb.A*[A_mat_pd+dyn.Ew*wext_x dyn.B+dyn.Ew*wext_u]];
          if dyn.np > 0
            Xd_b = [Xd_b; 
                    Xb.b-Xb.A*(F_mat_pd+dyn.Ew*(wext_p*PV(ip,:)'+wext_w))];
          else
            Xd_b = [Xd_b; 
                    Xb.b-Xb.A*(F_mat_pd+dyn.Ew*wext_w)];
          end
        end
      else
        Xd_A = [Xd_A; 
                Xb.A*[A_mat_pd dyn.B]];
        Xd_b = [Xd_b; 
                Xb.b-Xb.A*F_mat_pd];
      end
    end
    pre_proj = Polyhedron('A', [Xd_A; dyn.XU.A], ...
                          'b', [Xd_b; dyn.XU.b]);
    
    proj{iter} = projection(pre_proj, 1:dyn.nx);
    proj{iter}.minHRep;
  end

  X0 = Polyhedron('H', cell2mat(cellfun(@(p) p.H, proj, 'UniformOutput', false)));
  X0.minHRep;
end
