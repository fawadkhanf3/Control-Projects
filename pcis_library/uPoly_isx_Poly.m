function [intersx_array] = uPoly_isx_Poly( polyU , poly )
	%Description
	%	Creates the intersection between a union of polyhedra and a given polyhedron.
	%

	%% Input Processing

	if ~isa(polyU,'PolyUnion')
		error('First input is not a PolyUnion.')
	end

	%% Algorithm

	intersx_array = [];
	for i_pU = 1:polyU.Num
		
		temp_intersect = intersect(polyU.Set(i_pU),poly);

		if ~temp_intersect.isEmptySet
			%outputUnion.add( polyU.Set(i_pU) )
			intersx_array = [ intersx_array temp_intersect ];
		end
	end

end