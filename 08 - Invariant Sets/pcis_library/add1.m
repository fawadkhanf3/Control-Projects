function union = add1(poly1, poly2)
	union = PolyUnion;
	if isa(poly1, 'PolyUnion')
		for i=1:poly1.Num
			union.add(poly1.Set(i));
		end
	else
		union.add(poly1);
	end
	if isa(poly2, 'PolyUnion')
		for i=1:poly2.Num
			union.add(poly2.Set(i));
		end
	else
		union.add(poly2);
	end
end
