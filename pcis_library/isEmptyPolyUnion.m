function bool = isEmptyPolyUnion(X)

    if isa(X,'PolyUnion')
        bool = true;
        for i = 1:X.Num
            if ~isEmptySet(X.Set(i))
                bool = false;
                return;
            end
        end
    else
        bool = isEmptySet(X);
    end
end