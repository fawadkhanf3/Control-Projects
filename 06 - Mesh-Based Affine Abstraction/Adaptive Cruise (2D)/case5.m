function sep = case5(hess,bounds)

fun = @(x) -norm([hess{2,1}(x(1),x(2)),hess{1,2}(x(1),x(2));hess{2,1}(x(1),x(2)),hess{2,2}(x(1),x(2))],2);
A = [];
b = [];
Aeq = [];
beq = [];
lb = [bounds(1,1),bounds(2,1)];
ub = [bounds(1,2),bounds(2,2)];
x0 = [0,0];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
sep = -fun(x);
end
