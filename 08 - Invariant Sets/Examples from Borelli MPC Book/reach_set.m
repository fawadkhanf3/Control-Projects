function Rch = reach_set(dyn,Goal,Container)

i = 0;
Ci = [Goal];
maxiter = 100;

while i < maxiter
    disp(['iteration ', num2str(i)])
    
    C2 = intersect(Container,predecessor_set(dyn,Ci(end)));    
    
    if isEmptySet(C2\Ci(end))
        break;
    end
    
    plot(C2, 'color', 'blue', 'alpha', 0.1);
    drawnow;
    
    Ci = [Ci;C2];
    i = i+1;
end

Rch = Ci;

% C.plot