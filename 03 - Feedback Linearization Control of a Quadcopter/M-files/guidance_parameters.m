function fp = guidance_parameters(fp)

syms t

fp.xref     = @(t) 1/2*cos(t/2);
fp.xrefdot  = matlabFunction(diff(fp.xref(t),t),'Vars',t);
fp.xref2dot = matlabFunction(diff(fp.xrefdot(t),t),'Vars',t);
fp.xref3dot = matlabFunction(diff(fp.xref2dot(t),t),'Vars',t);
fp.xref4dot = matlabFunction(diff(fp.xref3dot(t),t),'Vars',t);

fp.yref     = @(t) 1/2*sin(t/2);
fp.yrefdot  = matlabFunction(diff(fp.yref(t),t),'Vars',t);
fp.yref2dot = matlabFunction(diff(fp.yrefdot(t),t),'Vars',t);
fp.yref3dot = matlabFunction(diff(fp.yref2dot(t),t),'Vars',t);
fp.yref4dot = matlabFunction(diff(fp.yref3dot(t),t),'Vars',t);

fp.zref     = @(t) 1+t/10;
fp.zrefdot  = matlabFunction(diff(fp.zref(t),t),'Vars',t);
fp.zref2dot = matlabFunction(diff(fp.zrefdot(t),t),'Vars',t);
fp.zref3dot = matlabFunction(diff(fp.zref2dot(t),t),'Vars',t);
fp.zref4dot = matlabFunction(diff(fp.zref3dot(t),t),'Vars',t);

fp.psiref     = pi/3;
fp.psirefdot  = 0;
fp.psiref2dot = 0;

end