clear;close all;clc

tauP1 = 2;
tauP2 = 0.5;
tau_sen = 0.05;
Kproc = 4;

denGp = conv([tauP1 1],[tauP2 1]);
Gp = tf(Kproc,denGp);
H = tf(1,[tau_sen 1]);
GpH = Gp*H;

%% 

Kp = 0.48;
Ki = 0.27;

s = tf('s');

PI = Kp + Ki/s;

[A,B,C,D] = ssdata(PI);

A_Hanus = A - (B/D)*C;
B_Hanus = [0 B/D];
C_Hanus = C;
D_Hanus = [D 0];
