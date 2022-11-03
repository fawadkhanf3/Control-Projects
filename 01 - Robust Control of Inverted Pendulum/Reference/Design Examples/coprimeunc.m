function Ks=coprimeunc(Gs,gamrel);

[a,b,c,d]=ssdata(Gs);
S=eye(size(d'*d))+d'*d;
R=eye(size(d*d'))+d*d';
Sinv=inv(S); Rinv=inv(R);
A1=a-b*Sinv*d'*c; R1=S; B1=b; Q1=c'*Rinv*c;
[X,XAMP,G]=care(A1,B1,Q1,R1);
A2=A1'; Q2=b*Sinv*b'; B2=c'; R2=R;
[Z,ZAMP,G]=care(A2,B2,Q2,R2);

XZ=X*Z;
gammin=sqrt(1+max(eig(XZ)))
gam=gamrel*gammin
gam2=gam*gam;
gamconst=(1-gam2)*eye(size(XZ));
Lc=gamconst+XZ; Li=inv(Lc'); Fc=-Sinv*(d'*c+b'*X);
Ac=a+b*Fc+gam2*Li*Z*c'*(c+d*Fc);
Bc=gam2*Li*Z*c';
Cc=b'*X;
Dc=-d';
Ks=ss(Ac,Bc,Cc,Dc);

