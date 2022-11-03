function lambda = lipchitz(deriv,bounds)

sdpvar xo yo zo;

constr = [ bounds(1,1)<= xo <= bounds(1,2), bounds(2,1)<= yo <= bounds(2,2), bounds(3,1)<= zo <= bounds(3,2)];

options = sdpsettings('verbose',0);
optimize(constr, -deriv{1,1}(xo,yo,zo),options); 
B1 = deriv{1,1}(xo,yo,zo);

optimize(constr, -deriv{1,2}(xo,yo,zo),options);
B2 = deriv{1,2}(xo,yo,zo);

optimize(constr, -deriv{1,3}(xo,yo,zo),options);
B3 = deriv{1,3}(xo,yo,zo);

lambda = sqrt(B1^2+B2^2+B3^2);

end