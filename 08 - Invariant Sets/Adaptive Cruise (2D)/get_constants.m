function con = get_constants()

con.mass = 1370; % Mass [kg]
con.g = 9.82; % Gravitational Acceleration [m/s^2]

con.f0    = 3.8e-3*con.mass*con.g; % Resistance Coefficient
con.f1    = 2.6e-5*con.mass*con.g; % Resistance Coefficient
con.f2    = 0.4161; % Resistance Coefficient

con.v_bar  = 15;
con.v_min  = 0;
con.v_max  = 35;
con.v_L    = 12;
con.v_des  = 34;

con.h_min   = 0;
con.h_max   = 300;
con.h_safe  = 3;
con.tau_des = 1.3;
con.tau_min = 1;

con.u_min = -0.3*con.mass*con.g;
con.u_max = +0.2*con.mass*con.g;

con.f0_bar = con.f0-con.f2*con.v_bar^2;
con.f1_bar = con.f1+2*con.f2*con.v_bar;

con.u_min_bar = con.u_min;
con.u_max_bar = con.u_max-con.f2*max((con.v_max-con.v_bar)^2,(con.v_min-con.v_bar)^2);

% Sample Time
con.dT = 0.5;