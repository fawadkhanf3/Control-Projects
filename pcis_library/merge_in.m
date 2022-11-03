function newpoly = merge_in(poly1, poly2)
	% MERGE_IN: Merge neighboring polyhedra.
	% ======================================================
	%
	% SYNTAX
	% ------
	%	U = merge_in(poly1,poly2)
	%
	% DESCRIPTION
	% -----------
	%	Computes a convex polyhedron contained in the union of two 
	%   neighboring polyhedra, s.t. all inequalities of poly1 are 
	%   preserved except the common facet.
	%
	% INPUT
	% -----
	%	poly1, poly2	Polyhedra representing the union
	% 				 	Class: Polyhedron
	% 	

	[tr, index1, index2] = isNeighbor(poly1,poly2);

	if ~tr
		error('merge_in: polyhedra not neighbors')
	end

	newH = poly1.H(setdiff(1:size(poly1.H, 1), index1), :);
	
	for i = 1:size(poly2.H, 1)
		sol = poly1.extreme(poly2.A(i,:)); % maximize in direction of new facet
		ximax = sol.supp;
		if ximax < poly2.b(i)
			newH = [newH; poly2.H(i,:)];
		end
	end

	newpoly = Polyhedron('H', newH);
end