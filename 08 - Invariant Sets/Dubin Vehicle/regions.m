function x_y_theta = regions(con,planes)

theta_min = planes{5}(1,1);
theta_max = planes{5}(1,2);

A = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1];
b = [con.x_max;-con.x_min;con.y_max;-con.y_min;theta_max;-theta_min];

x_y_theta = Polyhedron(A,b);