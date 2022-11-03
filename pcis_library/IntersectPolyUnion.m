function [P12] = IntersectPolyUnion(P1,P2, MaxNum)
    if nargin < 2
        error('Not enough input arguments.');
    elseif nargin == 2
        MaxNum = inf;
    end
    
    if isa(P1,'Polyhedron')
        P1 = PolyUnion(P1);
    end
    if isa(P2,'Polyhedron')
        P2 = PolyUnion(P2);
    end
    
    P12 = PolyUnion;
    
    for i = 1:P2.Num
        for j = 1:P1.Num
            int_ji = P1.Set(j).intersect(P2.Set(i));
            if ~int_ji.isEmptySet
                int_ji.minHRep;
                P12.add(int_ji);
            end
            
%             if(P12.Num >= MaxNum)
%                 P12.reduce; % or P12.merge
%                 MaxNum = max(MaxNum,P12.Num*2);
%             end
        end
    end
    
%     P12.reduce;
%     % reduce the number of the Polytopes
%     if(P12.Num >= MaxNum)
%         P12.reduce; % or P12.merge
%     end
end

