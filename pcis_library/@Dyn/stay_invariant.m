function C = stay_invariant( varargin )
% Find controlled invariant set (CIS) C \subseteq S
% G is CIS to be expanded 
%
% Usage:
%   - C = stay_invariant(dyn, S, G, rho)
%   - C = stay_invariant(dyn, S, G, rho, debug_flag )
%
% Inputs:
%   debug_flag: - "0" will produce no debugging information.
%               - "1" will produce some textual debugging information
%               - "2" will produce some figure-based/plots of debugging information

%% Process Inputs

dyn = varargin{1};
S   = varargin{2};
G   = varargin{3};
rho = varargin{4};

if nargin >= 5
    debug_flag = varargin{5};
else
    debug_flag = 0;
end

%% Start Algorithm

C = PolyUnion(G);
sets2explore = PolyUnion(G); % new additions to C
preSets = []; 
is_converged = 0;
iter = 0;
while ~is_converged     % repeat until C doesn't grow (can also do it for a fixed # steps)
    for i = 1:sets2explore.Num % for each new set
        pre_c = pre(dyn, sets2explore.Set(i), rho); % find pre
        for s = 1:S.Num % take the intersection of this pre with each safe set
            p = intersect(pre_c, S.Set(s));
            if ~p.isEmptySet() % discard if empty
                preSets = [preSets p];
            end
        end
    end
    if isempty(preSets)
        1; % for debugging purposes
    end
    
    % Union of all pre (simplified by merging whenever possible)
    preSets = PolyUnion(preSets);
    preSets = preSets.merge;
    
    % check convergence (if all pre is contained in C)
    is_contained = zeros(preSets.Num,1);
    for i = 1:preSets.Num
        pi = preSets.Set(i);
        pis = pi - Polyhedron('A', [eye(pi.Dim); -eye(pi.Dim)], 'b', repmat(rho,2*pi.Dim,1));
        pr = [];
        for j = 1:C.Num
            rj = C.Set(j);
            pr = [pr intersect(pi,rj)];
        end
        pr = PolyUnion(pr);
        pr = pr.merge;
        if pr.Num == 1
            pr = pr.Set(1);
            is_contained(i) =  pr.contains(pis);
        end
    end

    is_converged = sum(is_contained) == preSets.Num;
    if is_converged
        1; % for debugging purposes
    end
    
    if debug_flag > 0
        
        %Check to see that the preSets.Set(i) does not intersect with any
        %objects in C PolyUnion
        disp(['C.Num = ' num2str(C.Num)])
        disp(['preSets.Num = ' num2str(preSets.Num)])
        for i_C = 1:C.Num
            for i_pS = 1:preSets.Num
                temp_intersect = intersect(C.Set(i_C),preSets.Set(i_pS));
                if temp_intersect.isEmptySet
                    disp(['Intersection of C.Set(' num2str(i_C) ') and preSets.Set(' num2str(i_pS) ') is empty!' ])
                else
                    if debug_flag > 1
                        figure;
                        hold on;
                        plot(C.Set(i_C).projection([2 3]));
                        plot(preSets.Set(i_pS).projection([2 3]),'color','blue')
                    end
                end
            end
        end

    end

%     %Get rid of the parts of preSets that don't overlap.
%     preSets_nix = [];
%     for i_pS = 1:preSets.Num
%         s0 = preSets.Set(i_pS); 
% 
%         i_intersx = [];
%         for i_C = 1:C.Num
%             %Detect the sets for which s0 intersects C
%             temp_intersect = intersect(C.Set(i_C),preSets.Set(i_pS));
%             if ~temp_intersect.isEmptySet
%                 i_intersx = [ i_intersx i_C];
%             end
% 
%         end
% 
%         %Find the sets that are not in either intersection peice.
%         temp_sm = [];
%         for i_C = i_intersx
%             temp_sm = [ temp_sm s0 \ C.Set(i_C) ];
%         end
% 
%         %Add this to the set of Non Intersecting PreSets (preSets_nix), if it is not empty
%         for i_s0 = 1:length(temp_s0)
%             if ~s0.isEmptySet
%                 preSets_nix = [preSets_nix s0]
%             end
%         end
%     end

    % add new found sets to C
    figure; hold on;
    for i = 1:preSets.Num
        plot(preSets.Set(i).projection([2 3]))
        C.add( preSets.Set(i) );
    end
    C = C.merge;
    % explore this preSets in the next iteration
    sets2explore = preSets;
    preSets = [];

    if debug_flag > 0
        iter = iter + 1;
        disp([ num2str(iter) ' Iterations Completed.']); % debugging purposes
        disp(' ')
    end
end
1; % debugging purposes
        
        