clear all;clc
Nominal_Aircraft_Plant_page_198;
G=ss(A,B,C,D);
G=-minreal(G);
w=logspace(-3,3,1e3);
tzero(G)
damp(G)
figure(1),impulse(G);grid
figure(2),step(G);grid
figure(3),margin(G);
%return
H=bodeoptions;
H.MagScale='log';H.MagUnits='abs';
H1=bodeoptions;
H1.MagScale='linear';H1.MagUnits='abs';

M=1.5; A=1e-3; wB=2;
Wp=tf([1/M wB],[1 wB*A]);
tau=0.04; r0=0.2; rinf=2;
Wm=tf([tau r0],[tau/rinf 1]);

%%
W1=55*tf([1 1],[1 0])*tf([1 4],[1 60]);
W1=50*tf([1 1],[1 0])*tf([1 4],[1 60]);

%W1=400*tf([1 1],[1 0])*tf([1 4],[1 60]);
%W1=30*tf([1 1],[1 0])*tf([1 2],[1 20]);
%W1=250*tf([1 1],[1 0])*tf([1 2],[1 60]);
figure(4),bodemag(W1,w,H);grid
Gs=G*W1;
figure(5),bode(G,'b',Gs,'r',w,H);grid;legend('Plant','Shaped Plant')
figure(6),margin(Gs);
%return
%%
Ks=coprimeunc(Gs,1.05);
Ks=-Ks;
K=W1*Ks;
size(K)
figure(7),margin(G*K)
%return
figure(8),bode(Gs,'b',G*K,'r');grid;legend('Shaped Plant','Shaped Plant with H-inf Controller')
%return
%%
% Analysis
L=G*K;      % Nominal loop transfer function
S=inv(1+L); % Nominal Sensitivity
T=1-S;      % Nominal complementary sensitivity
figure(9),subplot(2,1,1),bodemag(S*G,'b',S,'k',1/Wp,'r.',w,H);legend('S*G','S','1/Wp')
figure(9),subplot(2,1,2),bodemag(Wp*S,'b',Wp*S*G,'r',w,H1);grid;title('NP'),legend('Wp*S','Wp*S*G')
figure(10),subplot(2,1,1),bodemag(T,'b',1/Wm,'r.',w,H);legend('T','1/Wm')
figure(10),subplot(2,1,2),bodemag(Wm*T,w,H1);grid;title('RS - Wm*T')
WpSG=bode(Wp*S*G,w);WpSG=squeeze(WpSG);
WT=bode(Wm*T,w);WT=squeeze(WT);
WpS=bode(Wp*S,w);WpS=squeeze(WpS);
figure(11),subplot(2,1,1),semilogx(w,WpSG+WT);grid;title('RP - input dist')
figure(11),subplot(212),semilogx(w,WpS+WT);grid;title('RP - output dist')
norm([Wp*S*G;Wm*T],'inf')
norm([Wp*S;Wm*T],'inf')
%return
%%
figure(312),subplot(2,1,1),step(T);grid;title('ref step response')
figure(312),subplot(212),step(K*S);grid
figure(313),subplot(2,1,1),step(S*G);grid;title('input dist step response')
figure(313),subplot(212),step(-K*S*G);grid
figure(314),subplot(2,1,1),step(S);grid;title('output dist step response')
figure(314),subplot(212),step(-K*S);grid
figure(15),bodemag(Wm*T,H);grid

%return
%%
Delta=ultidyn('Delta',[1 1]);
Gp=G*(1+Wm*Delta);
Lp=Gp*K;
Sp=inv(1+Lp);
Tp=1-Sp;
[StabMarg,DeStabUnc,Reps]=robuststab(Tp);
StabMarg
%zpk(DeStabUnc.Delta)
[PerfMargin,PerfMarginUnc,Repp]=robustperf(Wp*Sp*Gp);
PerfMargin
% The above doesn't match with the blk diagram since the dist enters in
% between perturbation and nominal plant G whereas here the dist is before
% the perturbed plant
%zpk(PerfMarginUnc.Delta)

[PerfMargin,PerfMarginUnc,Repp]=robustperf(Wp*Sp);
PerfMargin
%return
figure(316),step(usample(Tp,50));grid;title('Ref Step response')
figure(317),step(usample(Sp*Gp,50));grid;title('Input Dist Step response')
figure(318),step(usample(Sp,50));grid;title('Output Dist Step response')
figure(319),bodemag(Wp*usample(Sp*Gp,50),w,H1);grid;title('RP - Input Dist')
figure(320),bodemag(Wp*usample(Sp,50),w,H1);grid;title('RP - Output Dist')
%axis([1e-3 1e3 .1 1])
