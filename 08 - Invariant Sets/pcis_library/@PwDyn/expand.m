function P = expand(varargin)
% Find controlled invariant set (CIS) C \subseteq S
% G is CIS to be expanded 
%
% Usage:
%   - P = expand(pwd0, S, G, rhoPre)
%   - P = expand(pwd0, S, G, rhoPre, 'plot_stuff')
%   - P = expand(pwd0, S, G, rhoPre, 'debug')
%   - P = expand(pwd0, S, G, rhoPre, 'plot_stuff','debug')
%
% Inputs:
%   plot_stuff: string, indicates that the system will plot figures of intermediate sets while operating.

%% Manage Inputs

pwd0 = varargin{1};
S = varargin{2};
G = varargin{3};
rhoPre = varargin{4};

%Assign defaults to the extra variables.
plot_stuff = 1;% false;
debug_flag = 1;%false;

for arg_ind = 5:nargin
    switch varargin{arg_ind}
    case 'plot_stuff'
        plot_stuff = true;
    case 'debug'
        debug_flag = true;
    otherwise
        error(['Unrecognized string input: ' varargin{arg_ind} ])
    end
end

%% Main Function

converged = zeros(S.Num,1);
P = cell(S.Num,1);
for i = 1:S.Num
    P{i} = G;
end

iter_num = 0;
counter = 0;
while sum(converged) < S.Num   
    
    for i = 1:S.Num
        if converged(i) == 1
            continue
        end
        if plot_stuff
            plot(P{i});%(end));
            drawnow;
        end
        
        pwd_pre = pwd0.pre(P{i}(end), rhoPre);
        
        if isempty(pwd_pre)
            return;
        end
            
        pwd_pre = pwd_pre.merge;
        P_next = PolyUnion(uPoly_isx_Poly(pwd_pre,S.Set(i)));
        P_next = P_next.merge;
        P{i} = [P{i} P_next];
        
        if isEmptySet(mldivide1(P{i}(end),P{i}(end-1)))
           converged(i) = 1;
           if sum(converged) == S.Num
               return
           end
        end

        %Update Loop Counter
        iter_num = iter_num + 1;
        if debug_flag
            disp(['Iteration #' num2str(iter_num) ' Complete.'])
            disp(' ')
        end
    end
    counter = counter+1;
end
        