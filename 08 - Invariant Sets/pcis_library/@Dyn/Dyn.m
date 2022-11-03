classdef Dyn
  % Dyn(A, F, B, XU, Ap, Fp, P, Ad, Fd, D, Ev, XPV_V, Ew, XPW_V): 
  % Discrete-time system of the form
  %
  %  x(k+1) = (A + ∑d_i Ad{i} + ∑p_i Ap{i}) x(k) + B u(k) + F ...
  %           + ∑d_i(k) Fd{i} + ∑p_i(k) Fp{i} + Ev v(k) + Ew w(k)
  %
  % u input such that (x(k),u(k)) ∈ XU
  % p(k) ∈ P measurable disturbance
  % d(k) ∈ D non-measurable disturbance
  % v(k) ∈ conv_i(XPV_V{i} [x(k); p(k); 1]) measurable state-dependent disturbance
  % w(k) ∈ conv_i(XPW_V{i} [x(k); p(k); 1]) non-measurable state-dependent disturbance
  %
  % INPUTS (all except A are optional and can be set later)
  %
  % A: (nx x nx) matrix
  % F: (nx x 1) matrix
  % B: (nx x nu) matrix
  % Ev: (nx x nv) matrix
  % Ew: (nx x nw) matrix
  %
  % Ap: cell of np matrices of size (nx x nx)
  % Fp: cell of np matrices of size (nx x 1)  
  % Ad: cell of nd matrices of size (nx x nx)
  % Fd: cell of nd matrices of size (nx x 1) 
  %
  % XU: (nx+nu)-dim Polyhedron
  % P:  np-dim Polyhedron
  % D:  nd-dim Pokyhedron
  %
  % XV_V: cell of nv matrices of size (1 x nx+np+1)  
  % XW_V: cell of nw matrices of size (1 x nx+np+1)  
  properties
    A;
    F;
    B;
    XU;
    Ap;
    Fp;
    P;
    Ad;
    Fd;
    D;
    Ev;
    XV_V;
    Ew;
    XW_V;
  end

  methods
    % Constructor
    function d = Dyn(A, F, B, XU, Ap, Fp, P, Ad, Fd, D, Ev, XV_V, Ew, XW_V)

      nx = size(A,2);

      d.A = A;

      if nargin < 2 || isempty(F)
        d.F = zeros(nx,1);
      else
        d.F = F;
      end

      if nargin < 3 || isempty(B)
        d.B = zeros(nx,0);
        d.XU = Polyhedron('H', [zeros(1,nx) 1]);
      else
        d.B = B;
        d.XU = XU;
      end

      if nargin < 5 || isempty(Ap)
        d.Ap = {};
        d.Fp = {};
        d.P = Polyhedron; 
      else
        d.Ap = Ap;
        d.Fp = Fp;
        d.P = P;
      end

      if nargin < 8 || isempty(Ad)
        d.Ad = {};
        d.Fd = {};
        d.D = Polyhedron; 
      else
        d.Ad = Ad;
        d.Fd = Fd;
        d.D = D;
      end

      if nargin < 11 || isempty(Ev)
        d.Ev = zeros(nx,0);
        d.XV_V = {};
      else
        d.Ev = Ev;
        d.XV_V = XV_V;
      end

      if nargin < 14 || isempty(Ew)
        d.Ew = zeros(nx,0);
        d.XW_V = {};
      else
        d.Ew = Ew;
        d.XW_V = XW_V;
      end
    end
    
    function check(d)
      % Checks
      assert(size(d.A, 1) == size(d.A,2))
      assert(size(d.F, 1) == size(d.A,1))
      assert(size(d.B, 1) == size(d.A,1))
      assert(size(d.Ew, 1) == size(d.A,1))
      assert(size(d.Ev, 1) == size(d.A,1))

      assert(d.XU.Dim == size(d.A,2) + size(d.B,2))

      assert(length(d.Ap) == d.P.Dim)
      assert(length(d.Fp) == d.P.Dim)
      assert(length(d.Ad) == d.D.Dim)
      assert(length(d.Fd) == d.D.Dim)

      for i=1:d.P.Dim
        assert(size(d.Ap{i}, 1) == size(d.A,1))
        assert(size(d.Ap{i}, 2) == size(d.A,1))
        assert(size(d.Fp{i}, 1) == size(d.A,1))
        assert(size(d.Fp{i}, 2) == 1)
      end

      for i=1:d.D.Dim
        assert(size(d.Ad{i}, 1) == size(d.A,1))
        assert(size(d.Ad{i}, 2) == size(d.A,1))
        assert(size(d.Fd{i}, 1) == size(d.A,1))
        assert(size(d.Fd{i}, 2) == 1)
      end

      for i=1:length(d.XV_V)
        assert(size(d.Ev,2) == size(d.XV_V{i}, 1))
        assert(size(d.A,2)+d.np+1 == size(d.XV_V{i}, 2))
      end

      for i=1:length(d.XW_V)
        assert(size(d.Ew,2) == size(d.XW_V{i}, 1))
        assert(size(d.A,2)+d.np+1+1 == size(d.XW_V{i}, 2))
      end
    end

    function n = nx(d)
      n = size(d.A,2);
    end
    function n = nu(d)
      n = size(d.B,2);
    end
    function n = nd(d)
      n = length(d.Ad);
    end
    function n = np(d)
      n = length(d.Ap);
    end
    function n = nw(d)
      n = size(d.Ew,2);
    end
    function n = nv(d)
      n = size(d.Ev,2);
    end
    function X0 = pre(d, X, rho)
      d.check();
      if nargin < 3
        rho = 0;
      end
      
      if isa(X, 'PolyUnion')	
		X0 = PolyUnion;
        for i=1:X.Num
          new_poly = pre_proj(d, X.Set(i), rho);
          X0.add(new_poly);
        end
      else
        X0 = pre_proj(d, X, rho);
      end
    end
  end
end

