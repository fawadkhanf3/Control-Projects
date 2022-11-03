function Inv = invariant_set(dyn,G)

i = 0;
Ci = G;
maxiter = 100;

while i < maxiter
    disp(['iteration ', num2str(i)])
    
    Pre_set = predecessor_set(dyn,Ci);
    Ci_1 = intersect(G, Pre_set);
    
    plot(Ci_1, 'alpha', 0.4,'color','green');
    drawnow;
    
%     if isEmptySet(mldivide(Ci, Ci_1))
%         break;
%     end
    if Ci_1<G
        break;
    end

    Ci = Ci_1;
    i = i+1;
end

Inv = Ci;

% C.plot