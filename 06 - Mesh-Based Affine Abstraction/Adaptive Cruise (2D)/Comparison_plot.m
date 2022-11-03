load planes.mat

f0 = 51;
f1 = 1.2567;
f2 = 0.4342;
mass = 1370;
g0 = 9.82;

% Bounds
v_min = 0;
v_max = 35;

Fw_min = -0.3*mass*g0;
Fw_max = +0.2*mass*g0;

func = @(v,Fw) (1/mass).*(Fw-f0-f1.*v-f2.*v.^2);
hyperplanes = planes;

num_grids = size(hyperplanes,1);
figure_number = 1;
plot_all(hyperplanes,func,figure_number);

hold on
 
plot_m = 50;

v_bar = 10;
f0 = 51;
f1 = 1.2567;
f2 = 0.4342;
f0_bar = f0-f2*(v_bar^2);
f1_bar = f1+2*f2*v_bar;

gamma = f2*max((v_max-v_bar)^2,(v_min-v_bar)^2);

Fw_min_bar = Fw_min;
Fw_max_bar = Fw_max - gamma;

v = linspace(v_min,v_max,plot_m);
Fw = linspace(Fw_min_bar,Fw_max_bar,plot_m);
[X,Y] = meshgrid(v,Fw);

vdot = 1/mass*(Y-f0_bar-(f1_bar)*X);


mesh(X,Y,vdot,'edgecolor','g ','linewidth',1.5);
mesh(X,Y,vdot+gamma,'edgecolor','y ','linewidth',1.5);
mesh(X,Y,vdot-gamma,'edgecolor','m ','linewidth',1.5)
    