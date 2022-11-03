clear all;clc
Nominal_Aircraft_Plant_page_198;
G=ss(A,B,C,D);
G=minreal(G);
w=logspace(-3,3,1e3);
tzero(G)
damp(G)
figure(1),impulse(G);grid
figure(2),step(G);grid
figure(3),margin(G);
%%
%return
H=bodeoptions;
H.MagScale='log';H.MagUnits='abs';
H1=bodeoptions;
H1.MagScale='linear';H1.MagUnits='abs';

M=1.5; A=1e-3; wB=2;
Wp=tf([1/M wB],[1 wB*A]);
figure(4),bodemag(1/Wp,H);grid
%return

tau=0.04; r0=0.2; rinf=2;
Wm=tf([tau r0],[tau/rinf 1]);
figure(5),bodemag(1/Wp,'b',Wm,'r',H);grid;axis([1e-1 1e3 0 10])
%return
%%
systemnames = 'G Wp Wm';
inputvar = '[d;u]';
outputvar = '[Wp;Wm;-G-d]';
input_to_G = '[u]';
input_to_Wp = '[G+d]';
input_to_Wm = '[G]';
P=sysic;
size(P)
%return
%%
%[K,CL,GAM,INFO] = hinfsyn(P,1,1,'display','on','method','lmi');
[K,CL,GAM,INFO] = hinfsyn(P,1,1,'display','on');
size(K)
L=G*K;  % Nominal loop transfer function
S=inv(1+L); % Nominal Sensitivity
T=1-S;      % Nominal Complementary Sensitivity
figure(6),subplot(2,1,1),bodemag(S,'b',1/Wp,'r',H);
figure(6),subplot(2,1,2),bodemag(Wp*S,H1);title('Nominal Performance')
figure(7),subplot(2,1,1),bodemag(T,'b',1/Wm,'r',H);
figure(7),subplot(2,1,2),bodemag(Wm*T,H1);title('Robust Stability')
WpSG=bode(Wp*S*G,w);WpSG=squeeze(WpSG);
WT=bode(Wm*T,w);WT=squeeze(WT);
WpS=bode(Wp*S,w);WpS=squeeze(WpS);
figure(8),subplot(2,1,1),semilogx(w,WpSG+WT);title('RP - mult unc & input dist rejection')
figure(8),subplot(2,1,2),semilogx(w,WpS+WT);title('RP - mult unc & output dist rejection')
norm([Wp*S*G;Wm*T],'inf')
norm([Wp*S;Wm*T],'inf')
%return
%%
figure(109),subplot(2,1,1),step(T);grid;title('ref step response')
figure(109),subplot(2,1,2),step(K*S);grid
figure(110),subplot(2,1,1),step(S*G);grid;title('input dist step response')
figure(110),subplot(2,1,2),step(-K*S*G);grid
figure(111),subplot(2,1,1),step(S);grid;title('output dist step response')
figure(111),subplot(2,1,2),step(-K*S);grid
figure(12),bodemag(Wm*T,H1);grid
%%
Delta=ultidyn('Delta',[1 1]);
Gp=G*(1+Wm*Delta);  % Perturbed Plant
Lp=Gp*K;            % Perturbed loop T/F
Sp=inv(1+Lp);       % Perturbed Sensitivity
Tp=1-Sp;            % Perturbed Complementary Sensitivity
[StabMarg,DeStabUnc,Reps]=robuststab(Tp);
StabMarg
zpk(DeStabUnc.Delta)
norm(DeStabUnc.Delta,'inf')
%return
[PerfMargin,PerfMarginUnc,Repp]=robustperf(Wp*Gp*Sp);
PerfMargin
%zpk(PerfMarginUnc.Delta)
%norm(PerfMarginUnc.Delta,'inf')
%return
[PerfMargin,PerfMarginUnc,Repp]=robustperf(Wp*Sp);
PerfMargin
figure(113),step(usample(Tp,50));grid;title('Ref Step response')
figure(114),step(usample(Sp*Gp,50));grid;title('Input Dist Step response')
figure(115),step(usample(Sp,50));grid;title('Output Dist Step response')
figure(116),bodemag(Wp*usample(Sp*Gp,50),w,H);grid;title('RP - Input Dist')
figure(117),bodemag(Wp*usample(Sp,50),w,H);grid;title('RP - Output Dist')
%axis([1e-3 1e3 .1 1])
%return
%%
%==============================
% mix sen S/T/KS
%Wu=0.05*tf([1],[1]);
%Wu=.5*tf([1 0.0001],[1 30]);
%%
Wu=.05*tf([1 0.1],[1 100]);
figure(200),bodemag(1/Wu,w,H);
%%
systemnames = 'G Wp Wm Wu';
inputvar = '[d;u]';
outputvar = '[Wp;Wm;Wu;-G-d]';
input_to_G = '[u]';
input_to_Wp = '[G+d]';
input_to_Wm = '[G]';
input_to_Wu = '[u]';
P=sysic;
size(P)
[K,CL,GAM,INFO] = hinfsyn(P,1,1,'display','on');
%[K,CL,GAM,INFO] = mixsyn(G,Wp,Wu,W);
%break
size(K)
L=G*K;  % loop transfer function
S=inv(1+L); % Sensitivity
T=1-S;      % complementary sensitivity
KS=K*S;
figure(26),subplot(2,1,1),bodemag(S,'b',1/Wp,'r',H);
figure(26),subplot(2,1,2),bodemag(Wp*S,H1);
figure(27),subplot(2,1,1),bodemag(T,'b',1/Wm,'r',H);
figure(27),subplot(2,1,2),bodemag(Wm*T,H1);
figure(28),subplot(2,1,1),bodemag(KS,'b',1/Wu,'r',H);
figure(28),subplot(2,1,2),bodemag(Wu*KS,H1);
WpSG=bode(Wp*S*G,w);WpSG=squeeze(WpSG);
WT=bode(Wm*T,w);WT=squeeze(WT);
WpS=bode(Wp*S,w);WpS=squeeze(WpS);
figure(29),subplot(2,1,1),semilogx(w,WpSG+WT);title('RP - mult unc & input dist rejection')
figure(29),subplot(2,1,2),semilogx(w,WpS+WT);title('RP - mult unc & output dist rejection')
norm([Wp*S*G;Wm*T;Wu*KS],'inf')
norm([Wp*S;Wm*T;Wu*KS],'inf')
%return

figure(209),subplot(2,1,1),step(T);grid;title('ref step response')
figure(209),subplot(2,1,2),step(K*S);grid
figure(210),subplot(2,1,1),step(S*G);grid;title('input dist step response')
figure(210),subplot(2,1,2),step(-K*S*G);grid
figure(211),subplot(2,1,1),step(S);grid;title('output dist step response')
figure(211),subplot(2,1,2),step(-K*S);grid

figure(30),bodemag(Wm*T,H1);grid
