function sep = case5(hess,bound)

fun = @(x) -norm(hess{1,1}(x),2);
A = [];
b = [];
Aeq = [];
beq = [];
lb = bound(1,1);
ub = bound(1,2);
x0 = [0,0];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
sep = -fun(x);
end
