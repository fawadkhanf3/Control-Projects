function mesh =  sample(m,bounds)

% [X,Y]=meshgrid(linspace(bounds(1,1),bounds(1,2),m),linspace(bounds(2,1),bounds(2,2),m));
% mesh=[X(:), Y(:), ones(m*m,1)];

Y = meshgrid(linspace(bounds(1,1),bounds(1,2),m));

mesh = [Y(:), ones(m*m,1)];
end