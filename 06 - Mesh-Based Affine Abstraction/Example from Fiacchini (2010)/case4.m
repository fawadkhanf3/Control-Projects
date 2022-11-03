function sep = case4(deriv,bounds)
sdpvar xo yo zo;

constr = [ bounds(1,1)<= xo <= bounds(1,2), bounds(2,1)<= yo <= bounds(2,2), bounds(3,1)<= zo <= bounds(3,2)];

options = sdpsettings('verbose',0);
obj=((deriv{1,1}(xo,yo,zo))^2+(deriv{1,2}(xo,yo,zo))^2+(deriv{1,3}(xo,yo,zo))^2);
optimize(constr, -obj,options); 

sep = sqrt(value(obj));

end