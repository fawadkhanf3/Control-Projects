function [xref,yref,zref,psiref] = reference_trajectory(t,fp)

xref(1) = fp.xref(t);
xref(2) = fp.xrefdot(t);
xref(3) = fp.xref2dot(t);
xref(4) = fp.xref3dot(t);

yref(1) = fp.yref(t);
yref(2) = fp.yrefdot(t);
yref(3) = fp.yref2dot(t);
yref(4) = fp.yref3dot(t);

zref(1) = fp.zref(t);
zref(2) = fp.zrefdot(t);
zref(3) = fp.zref2dot(t);
zref(4) = fp.zref3dot(t);

psiref(1) = fp.psiref;
psiref(2) = fp.psirefdot;

end