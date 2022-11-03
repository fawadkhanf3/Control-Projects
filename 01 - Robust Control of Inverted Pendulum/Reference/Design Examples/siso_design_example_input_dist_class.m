
Nominal_Aircraft_Plant_page_198;
G=ss(A,B,C,D);
G=minreal(G);
w=logspace(-3,3,1e3);
tzero(G)
damp(G)
figure(1),impulse(G);grid
figure(2),step(G);grid
figure(3),margin(G);

H = bodeoptions;
H.MagScale = 'log'; 
H.MagUnits = 'abs';

H1 = bodeoptions;
H1.MagScale = 'linear';
H1.MagUnits ='abs';

M  = 1.5; A=1e-3; wB=2;
Wp = tf([1/M wB],[1 wB*A]);

figure(4),bodemag(1/Wp,H);

tau = 0.04; r0=0.2; rinf=2;
Wm=tf([tau r0],[tau/rinf 1]);
figure(5),bodemag(1/Wp,'b',Wm,'r',H);grid;axis([1e-1 1e3 0 10])
%return

systemnames = 'G Wp Wm';
inputvar    = '[d;u]';
outputvar   = '[Wp;Wm;-G]';
input_to_G  = '[u+d]';
input_to_Wp = '[G]';
input_to_Wm = '[u]';
P=sysic;
size(P)
%%
%[K,CL,GAM,INFO] = hinfsyn(P,1,1,'display','on','method','lmi');
[K,CL,GAM,INFO] = hinfsyn(P,1,1,'display','on');
size(K)
L=G*K;  % Nominal loop transfer function
S=inv(1+L); % Nominal Sensitivity
T=1-S;      % Nominal complementary sensitivity
figure(6),subplot(2,1,1),bodemag(S*G,'b',1/Wp,'r',H);     % NP
figure(6),subplot(2,1,2),bodemag(Wp*S*G,H1);title('NP')
figure(7),subplot(2,1,1),bodemag(T,'b',1/Wm,'r',H);       % RS
figure(7),subplot(2,1,2),bodemag(Wm*T,H1);title('RS')

WpSG=bode(Wp*S*G,w);WpSG=squeeze(WpSG);
WT=bode(Wm*T,w);WT=squeeze(WT);
WpS=bode(Wp*S,w);WpS=squeeze(WpS);
figure(8),subplot(211),semilogx(w,WpSG+WT);title('RP - mult unc & input dist rejection')
figure(8),subplot(212),semilogx(w,WpS+WT);title('RP - mult unc & output dist rejection')
norm([Wp*S*G;Wm*T],'inf')
norm([Wp*S;Wm*T],'inf')

%%
figure(209),subplot(211),step(T);grid;title('ref step response')
figure(209),subplot(212),step(K*S);grid
figure(210),subplot(211),step(S*G);grid;title('input dist step response')
figure(210),subplot(212),step(-K*S*G);grid
figure(211),subplot(2,1,1),step(S);grid;title('output dist step response')
figure(211),subplot(212),step(K*S);grid

%%
figure(12),bodemag(Wm*T,H);grid
Delta=ultidyn('Delta',[1 1]);
Gp=G*(1+Wm*Delta);  % perturbed plant
Lp=Gp*K;            % pertured loop T/F
Sp=inv(1+Lp);       % perturbed Sensitivity
Tp=1-Sp;
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
figure(213),step(usample(Tp,50));grid;title('Ref Step response')
figure(214),step(usample(Sp*Gp,50));grid;title('Input Dist Step response')
figure(215),step(usample(Sp,50));grid;title('Output Dist Step response')
figure(216),bodemag(Wp*usample(Sp*Gp,50),w,H);grid;title('RP - Input Dist')
figure(217),bodemag(Wp*usample(Sp,50),w,H);grid;title('RP - Output Dist')
%axis([1e-3 1e3 .1 1])
