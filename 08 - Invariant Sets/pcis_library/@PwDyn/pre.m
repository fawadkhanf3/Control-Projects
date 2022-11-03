function [preSets] = pre( varargin )

	% Usage:
	%	preSets = pre(pwd1 , X , rho)
	%	preSets = pre(pwd1 , X , rho , 'plot_stuff')
	% Inputs:
	%	pwd1	- The piecewise affine dynamics we are using.
	%	X 		- Target Set. The set for which we are calculating the pre.
	%
	% Outputs:
	%	preSets - An array of polyhedra. All possible polyhedra that could result in the target states X.

	%% Input Processing

	pwd1 = varargin{1};
	X = varargin{2};
	rho = varargin{3};

	plot_flag = 0;

	if nargin >= 4
        arg_idx = 4;
        while arg_idx <= nargin
			switch varargin{arg_idx}
				case 'plot_stuff'
                    if nargin == arg_idx
                        plot_flag = 0;
                        arg_idx = arg_idx + 1;
                    elseif isnumeric(varargin{arg_idx+1})
                        plot_flag = varargin{arg_idx+1};
                        arg_idx = arg_idx + 2;
                    else
                        plot_flag = 1;
                        arg_idx = arg_idx + 1;
                    end
                    figure;
				otherwise
					error(['Unrecognized input: ' varargin{ind} ])
			end
		end
	end

    if isa(X,'PolyUnion')
        %If the target set is a polyUnion, then we need to do the pre for
        %each element in its 'Sets'.
        temp_preSets = cell(X.Num,1);
        for ind_X = 1:X.Num
            temp_preSets{ind_X} = pwd1.pre(X.Set(ind_X),rho);
        end
        %After all is said and done, combine the polyUnions.
        all_sets = [];
        for ind_X = 1:length(temp_preSets)
            all_sets = [all_sets temp_preSets{ind_X}.Set ];
        end
        
        if ~isempty(all_sets)
            preSets = PolyUnion(all_sets);
        else
            warning('Pre contained only empty sets?')
            preSets = [];
        end
        
    elseif isa(X,'Polyhedron')

	% eps0 = 0.01;

	preSets = PolyUnion;
	for i=1:pwd1.num_region
	    new_poly = pwd1.reg_list{i}.intersect(pwd1.dyn_list{i}.pre_proj(X, rho));
	    % S = add1(S, new_poly);		
		preSets.add(new_poly);
	end

end