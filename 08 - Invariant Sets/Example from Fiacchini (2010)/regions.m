function x_y_theta = regions(con,planes) %#ok<*INUSL>

x1_min = planes{5}(1,1);
x1_max = planes{5}(1,2);

x2_min = planes{5}(2,1);
x2_max = planes{5}(2,2);

A = [1 0;-1 0;0 1;0 -1];
b = [x1_max;-x1_min;x2_max;-x2_min];

x_y_theta = Polyhedron(A,b);