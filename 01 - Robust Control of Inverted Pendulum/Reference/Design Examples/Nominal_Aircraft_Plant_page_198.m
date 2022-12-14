clc
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Nominal Aircraft Plant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A=[0 0 1.1320 0 -1;
    0 -0.0538 -0.1712 0 0.0705;
    0 0 0 1 0;
    0 0.0485 0 -0.8556 -1.0130;
    0 -0.2909 0 1.0532 -0.6859];
B=[0 0 0;-0.12 1 0;0 0 0;4.4190 0 -1.6650;1.5750 0 -0.0732];
C=[eye(3) zeros(3,2)];
D=zeros(3);
B=B(:,3);
C=C(3,:);
D=0;

%control input=elevator angle (deg)
%output=pitch angle (deg)