function sep = case4(deriv,bounds)
sdpvar yo;

constr = [bounds(1,1)<= yo <= bounds(1,2)];

options = sdpsettings('verbose',0);
obj= (deriv{1,1}(yo))^2;
optimize(constr, -obj,options); 

sep = sqrt(value(obj));

end