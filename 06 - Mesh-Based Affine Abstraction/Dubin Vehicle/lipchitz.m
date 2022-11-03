function lambda = lipchitz(deriv,bounds)

sdpvar xo;

constr = bounds(1,1)<= xo <= bounds(1,2);

options = sdpsettings('verbose',0);
optimize(constr, -deriv{1,1}(xo),options); 
lambda = deriv{1,1}(xo);

end