function lambda = lipchitz(deriv,bounds)

sdpvar xo yo;

constr = [ bounds(1,1)<= xo <= bounds(1,2), bounds(2,1)<= yo <= bounds(2,2)];

options = sdpsettings('verbose',0);
optimize(constr, -deriv{1,1}(xo,yo),options); 
B1 = deriv{1,1}(xo,yo);

optimize(constr, -deriv{1,2}(xo,yo),options);
B2 = deriv{1,2}(xo,yo);

lambda = sqrt(B1^2+B2^2);

end