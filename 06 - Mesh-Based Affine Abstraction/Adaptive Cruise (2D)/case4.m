function sep = case4(deriv,bounds)
sdpvar xo yo;

constr = [ bounds(1,1)<= xo <= bounds(1,2), bounds(2,1)<= yo <= bounds(2,2)];

options = sdpsettings('verbose',0);
obj=((deriv{1,1}(xo,yo))^2+(deriv{1,2}(xo,yo))^2);
optimize(constr, -obj,options); 

sep = sqrt(value(obj));

end