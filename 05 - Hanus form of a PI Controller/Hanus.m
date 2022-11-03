%% Hanus form for PI Controller
clear;close all;
format long g

Kp = 8e-5;
Ki = Kp/10;

s = tf('s');

PI = Kp + Ki/(s+0.01);

[A,B,C,D] = ssdata(PI);

A_Hanus = A - (B/D)*C;
B_Hanus = [0 B/D];
C_Hanus = C;
D_Hanus = [D 0];

sys_Hanus  = ss(A_Hanus,B_Hanus,C_Hanus,D_Hanus);
sys_HanusD = c2d(sys_Hanus,0.02,'Tustin');

[A,B,C,D] = ssdata(sys_HanusD);
