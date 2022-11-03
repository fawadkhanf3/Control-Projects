function Pre = predecessor_set(dyn,G)

n = dyn.n;

Ad = dyn.A;
Bd = dyn.B;
Ed = dyn.E;
Kd = dyn.K;

Hx = G.A;
hx = G.b;
Hu = dyn.Hu; % [Hux Huu];
hu = dyn.hu;

Hdx_minus = dyn.Hdx_minus;
Hdu_minus = dyn.Hdu_minus;

Hdx_plus = dyn.Hdx_plus;
Hdu_plus = dyn.Hdu_plus;

HH = [Hx*(Ad+Ed*Hdx_plus)  Hx*Bd;
      Hx*(Ad+Ed*Hdx_minus) Hx*Bd;
      Hu];

hh = [hx-Hx*(Ed*Hdu_plus+Kd);
      hx-Hx*(Ed*Hdu_minus+Kd);
      hu];

P = Polyhedron(HH,hh);
Pre = P.projection(1:n);