function mesh =  sample(m,bounds)

[X,Y,Z] = meshgrid(linspace(bounds(1,1),bounds(1,2),m),linspace(bounds(2,1),bounds(2,2),m),linspace(bounds(3,1),bounds(3,2),m));

mesh = [X(:), Y(:), Z(:), ones(m*m*m,1)];

end