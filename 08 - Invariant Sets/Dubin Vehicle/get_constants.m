function con = get_constants()

con.theta_min = 0;%-0.44;
con.theta_max = 2*pi;%0.44;

con.x_min = 0;
con.x_max = 200;
con.y_min = 0;
con.y_max = 200;

con.u_min = 0;%-tan(con.theta_max);%-0.5*pi/180;
con.u_max = 2*pi/180;%tan(con.theta_max);%v/r*tand(2);%+0.5*pi/180;

con.dT = 0.5;

