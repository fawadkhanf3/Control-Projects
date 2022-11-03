function [X0] = pre_rc(dyn, X, rho)
  % Compute inner approximation of set
  % {x : ∀ p ∃ u ∀ d, x(t+1) + Ball(rho) ⊆ X}
  % using robust counterpart method
  %
  % Reference: Petter Nilsson Ph.D. thesis (2017), Theorem 3.5

  if ~isa(dyn, 'Dyn')
    error('dyn must be an instance of Dyn');
  end

  if nargin < 3
    rho = 0;
  end

  if length(rho) == 1
    rho = rho * ones(dyn.nx,1);
  end

  Xb = X - Polyhedron('A', [eye(dyn.nx); -eye(dyn.nx)], 'b', repmat(rho,2,1));
  Xb.minHRep;   % important

  if norm(dyn.XU.A(:, 1:dyn.nx)) > 0
    error('state-dependent disturbance not supported yet')
  end

  if norm(dyn.Ew) > 0
    error('E w disturbance not supported yet')
  end

  if norm(dyn.Ev) > 0
    error('E w disturbance not supported yet')
  end

  X0 = pre_mat(dyn.A, dyn.Ap', dyn.Ad', dyn.B, dyn.Fp', dyn.Fd', ...
           Xb.A, Xb.b, dyn.P.A, dyn.P.b, dyn.D.A, dyn.D.b, ...
           dyn.XU.A(:, dyn.nx+1:end), dyn.XU.b);
end


function [proj] = pre_mat(A, A_d1, A_d2, B, F_d1, F_d2, ...
              H_x, h_x, H_d1, h_d1, ...
              H_d2, h_d2, H_u, h_u)
% pre: construct pre matrices for parametrized systems
%      in particular, construct Heq, heq, Hiq, hiq such that
% 
%      proj_x {Heq [x;y] = heq, Hieq [x;y] \leq hiq}
%      is contained in the set
% 
%     { x : \forall d1 \exists u \forall d2, H_x x^+ \leq x },
%
%      where
%
%    x^+ = (A + \sum_j d1_j A_d1^j + \sum_j d2_j A_d2^j)x + 
%        Bu + \sum_j d1_j F_d1^j + \sum_j d2_j F_d2^j
%
%      and for constraints  
%
%    H_u u \leq u, H_d1 d1 \leq h_d1, H_d2 d2 \leq h_d2
%
  [N_x, n_x] = size(H_x);
  [N_u, n_u] = size(H_u);
  [N_d1, n_d1] = size(H_d1);
  [N_d2, n_d2] = size(H_d2);

  % Some sanity checks
  assert(length(A_d1) == n_d1);
  assert(length(F_d1) == n_d1);
  assert(length(A_d2) == n_d2);
  assert(length(F_d2) == n_d2);

  n1 = n_x;         % x 
  n2 = N_d2 * N_x;      % vec(y0)
  n3 = N_d2 * N_x * n_d1;   % vec([Ky_1; ...; Ky_{n_d1}])
  n4 = n_u;         % vec(Ku)
  n5 = n_u * n_d1;      % vec(z1)
  n6 = N_d1 * N_x;
  n7 = N_d1 * N_x * n_d2;   % [vec(z2_1); ...; vec(z2_nd2)]
  n8 = N_d1 * N_x * n_d2;   % [vec(z3_1); ...; vec(z3_nd2)]
  n9 = N_d1 * N_x * N_d2;   % [vec(z4_1); ...; vec(z4_Nd2)]
  n10 = N_d1 * N_u;     % vec(z_5)

  n_tot = n1 + n2 + n3 + n4 + n5 + n6 + ...
        n6 + n7 + n8 + n9 + n10;

  % first ineq (N_x)
  m1 = N_x;
  H1_1 = H_x * A;
  H1_2 = kron(eye(N_x), h_d2');
  H1_4 = H_x * B;
  H1_6 = kron(eye(N_x), h_d1');

  h1 = vec(h_x);
  H1 = [H1_1 H1_2 zeros(m1, n3) H1_4 zeros(m1, n5) ...
      H1_6 zeros(m1, n7+n8+n9+n10)];

  % second eq (n_d1 * N_x)
  m2 = n_d1 * N_x;
  H2_1 = -vec_pem(N_x, n_d1) * kron(eye(n_d1), H_x) * cell2mat(A_d1);
  H2_3 = -kron(eye(N_x), kron(eye(n_d1), h_d2'));
  H2_5 = -vec_pem(N_x, n_d1) * kron(eye(n_d1), H_x*B);
  H2_6 = kron(eye(N_x), H_d1');

  if isempty(H2_1)
    H2_1 = zeros(m2, n1);
  end

  h2 = vec((kron(eye(n_d1), H_x) * cell2mat(F_d1))');
  H2 = [H2_1 zeros(m2, n2), H2_3 zeros(m2, n4) H2_5 H2_6 ...
      zeros(m2, n7+n8+n9+n10)];


  % third ineq (N_x * n_d2)
  m3 = N_x * n_d2;
  H3_1_list = {};

  H3_2_list = {};
  H3_7_list = {};
  h3_list = {};

  for i=1:n_d2
    H3_1_list{end+1} = -vec_pem(N_x,1) * kron(unit_vec(i, n_d2)', eye(N_x)) ...
      * vec_pem(N_x, n_d2) * kron(eye(n_d2), H_x) ...
      * cell2mat(A_d2);
    H3_2_list{end+1} = kron(eye(N_x), unit_vec(i, n_d2)'*H_d2');
    H3_7_list{end+1} = kron(eye(N_x), h_d1');

    h3_list{end+1} = vec_pem(N_x,1) * kron(unit_vec(i, n_d2)', eye(N_x)) ...
      * vec_pem(N_x, n_d2) * kron(eye(n_d2), H_x) ...
      * cell2mat(F_d2); 
  end

  H3_1 = cell2mat(H3_1_list');
  H3_2 = cell2mat(H3_2_list');
  H3_7 = blkdiag(H3_7_list{:});

  if isempty(H3_1)
    H3_1 = zeros(m3, n1);
  end

  h3 = cell2mat(h3_list');
  H3 = [H3_1 H3_2 zeros(m3, n3+n4+n5+n6) H3_7 zeros(m3, n8+n9+n10)];

  % fourth eq (n_d1 * N_x * n_d2)
  m4 = n_d1 * N_x * n_d2;
  H4_3_list = {};
  H4_7_list = {};
  h4_list = {};
  for i=1:n_d2
    H4_3_list{end+1} = -kron(eye(N_x), kron(eye(n_d1), unit_vec(i, n_d2)'*H_d2'));
    H4_7_list{end+1} = kron(eye(N_x), H_d1');
    h4_list{end+1} = zeros(n_d1*N_x, 1);
  end

  H4_3 = cell2mat(H4_3_list');
  H4_7 = blkdiag(H4_7_list{:});

  h4 = cell2mat(h4_list');
  H4 = [zeros(m4, n1+n2) H4_3 zeros(m4, n4+n5+n6) H4_7 zeros(m4, n8+n9+n10)];

  % fifth ineq (N_x * n_d2)
  m5 = m3;
  H5_1 = -H3_1;
  H5_2 = -H3_2;
  H5_8 = H3_7;

  h5 = -h3;
  H5 = [H5_1 H5_2 zeros(m3, n3+n4+n5+n6+n7) H5_8 zeros(m3, n9+n10)];

  % sixth eq (n_d1 * N_x * n_d2)
  m6 = m4;
  H6_3 = -H4_3;
  H6_8 = H4_7;

  h6 = h4;
  H6 = [zeros(m6, n1+n2) H6_3 zeros(m6, n4+n5+n6+n7) H6_8 zeros(m6, n9+n10)];

  % 7th ineq (N_x * N_d2)
  m7 = N_x*N_d2;
  H7_2_list = {};
  H7_9_list = {};
  h7_list = {};

  for i=1:N_d2
    H7_2_list{end+1} = -kron(eye(N_x), unit_vec(i, N_d2)');
    H7_9_list{end+1} = kron(eye(N_x), h_d1');
    h7_list{end+1} = zeros(N_x, 1);
  end

  H7_2 = cell2mat(H7_2_list');
  H7_9 = blkdiag(H7_9_list{:});

  h7 = cell2mat(h7_list');
  H7 = [zeros(m7, n1) H7_2 zeros(m7, n3+n4+n5+n6+n7+n8) H7_9 zeros(m7, n10)];

  % 8th eq (n_d1 * N_x * N_d2)
  m8 = n_d1 * N_x * N_d2;
  H8_3_list = {};
  H8_9_list = {};
  h8_list = {};

  for i=1:N_d2
    H8_3_list{end+1} = kron(eye(N_x), kron(eye(n_d1), unit_vec(i, N_d2)'));
    H8_9_list{end+1} = kron(eye(N_x), H_d1');
    h8_list{end+1} = zeros(n_d1 * N_x, 1);
  end

  H8_3 = cell2mat(H8_3_list');
  H8_9 = blkdiag(H8_9_list{:});

  h8 = cell2mat(h8_list');
  H8 = [zeros(m8, n1+n2) H8_3 zeros(m8, n4+n5+n6+n7+n8) H8_9 zeros(m8, n10)];

  % 9th ineq (N_u)
  m9 = N_u;
  H9_4 = H_u;
  H9_10 = kron(eye(N_u), h_d1');

  H9 = [zeros(m9, n1+n2+n3) H9_4 zeros(m9, n5+n6+n7+n8+n9) H9_10];
  h9 = vec(h_u');

  % 10th eq (n_d1 * N_u)
  m10 = n_d1 * N_u;
  H10_5 = -vec_pem(N_u, n_d1) * kron(eye(n_d1), H_u);
  H10_10 = kron(eye(N_u), H_d1');

  H10 = [zeros(m10,n1+n2+n3+n4) H10_5 zeros(m10,n6+n7+n8+n9) H10_10];
  h10 = zeros(m10, 1);

  % Positivity (N_d1 * N_x + 2*N_d1*N_x*n_d2 + N_d1*N_x*N_d2 + N_d1*N_u)
  m11 = N_d1 * N_x + 2*N_d1*N_x*n_d2 + N_d1*N_x*N_d2 + N_d1*N_u;
  H11_6_10 = -eye(m11);

  h11 = zeros(m11, 1);
  H11 = [zeros(m11, n1+n2+n3+n4+n5) H11_6_10];

  Heq = [H2; H4; H6; H8; H10];
  heq = [h2; h4; h6; h8; h10];

  Hiq = [H1; H3; H5; H7; H9; H11];
  hiq = [h1; h3; h5; h7; h9; h11];

  P = Polyhedron('H', [Hiq hiq], 'He', [Heq heq]);
  P.minHRep;

  proj = projection(P, 1:n_x);
end

%% vec_pem: compute vector permutation P such that
% vec(A') = P vec(A) for a matrix A of size (m x n) 
function [P] = vec_pem(m, n)
  P = zeros(m*n, m*n);

  for i=1:m
    ei = zeros(1,m);
    ei(i) = 1;
    P(1+(i-1)*n:i*n, :) = kron(eye(n), ei);
  end
end

%% unit_vec: n-dimensional i:th column unit vector
function [ei] = unit_vec(i, n)
  ei = zeros(n,1);
  ei(i) = 1;
end
