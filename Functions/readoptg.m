function [tArgs, ntArgs, Args] = readoptg(args, varargin)

% Read varargins as a list of pairs (fieldname, value)
% and creates a structure out of it.
% The first parameter, "arg", is a structure whose field names
% are the expected field names.
%
% - The first output "tArgs" is a structure that contains exactly
%   the same field names as args. For each field name, if it
%   is present in varargin, then the associated value is the one passed
%   in varargin, otherwise, the value is that stored in args.
% - The second output is a structure composed of the field that
%   where passed in varargin but were not present in args.
% - The third argument is the merged structure
%
% To increase the handeness, some additional feature are enabled which
% form an exception to the (x = fieldname, v = value) pair rule. These are :
%
% if
%     i)   x starts with ?, or
%     ii)  x starts with !, or
%     iii) x contains a white space or a equal sign
% then y is not passed (the next argument is not of "type" value but of
% type "fieldname") but it is insteasd encoded in x directly.
% In case i)-ii), if x starts with ? (resp. !), y takes the value (boolean)
% true (resp. false). y is therefore a boolean and cannot be something
% else.
% In case iii), then the field name is the fragment that preceeds : and y
% the fragment that follows. y is therefor a string (and cannot be
% something else).
%
% As an example the two calls are equivalent :
% pparams(arg, 'param1', true, 'param2', false, 'param3', 'abcd')
% pparams(arg, '?param1', '!param2', 'param3:abcd')
%
% Another feature is the '-' prefix. If x starts with a dash,
% the fieldname of arg that has the susequent keyword has a prefix if taken.
% An error is thrown if there is 0 or more than 1 match.
%
% Another feature is the case where the first optional argument is '#'. Is that case
% the list that follows is a list of values, the nth value is being associated with the nth field of
% the structure arg (the first argument passed to pparams)
%
% As an example, the two calls are equivalent :
%
% pparams(arg, 'abcd', 1, 'defg', true)
% pparams(arg, '-a', 1, '-?d')
%
% as long as there is no ambiguity with the shortcuts 'a' and 'b'
% (this may for ex. fail if arg contained a field named 'ab' in addition to 'abcd' and 'defg')

argf = fieldnames(args);

Args = struct();

if (numel(varargin) >= 1 && ischar(varargin{1}) && strcmp(varargin{1}, '#'))
    
    for k = 1:numel(varargin)-1
        
        if (k > numel(argf))
            break;
        end
        
        Args.(argf{k}) = varargin{k+1};
        
    end
    
elseif (numel(varargin) >= 1 && ischar(varargin{1}) && any(regexp(varargin{1}, '^#.+')))
    
    
    
    for k = 1:numel(varargin)
        
        if (k > numel(argf))
            break;
        end
        
        if (k == 1)
            Args.(argf{k}) = varargin{k}(2:end);
        else
            Args.(argf{k}) = varargin{k};
        end
        
        
    end
    
    
else
    
    % Fill-in Args, field by field (by looping through the varargin list)
    
    k = 1;
    nparam = 0;
    
    while (k <= numel(varargin))
        
        param = varargin{k};
        nparam = nparam + 1;
        
        % If param is a structure, nothing to do (already has the "output" format)
        
        if (isstruct(param))
            
            % Loop over the fields and checks wether field name already existed
            % If this occurs, warn and overwrite.
            
            f = fieldnames(param);
            
            for j = 1:numel(f)
                
                if (isfield(Args, f{j}))
                    
                    % parameter appears more than once
                    warning ('pTXObj.readopt', 'The argument ''%s'' appears more that once (the last passed value is retained!)', f{j});
                    
                end
                
                Args.(f{j}) = param.(f{j});
                
            end
            
            % If params is not a string, bad argument
            
        elseif (~ischar(param) || isempty(regexp(param, '[a-zA-Z]', 'once')))
            
            error('Parameter nb %d must be either a structure or character array containing a valid field name', k);
            
        else
            
            if (regexp(param, '^-'))
                
                % find in argf the string that starts with param
                
                param = findmatch(param, argf);
                
            end
            
            % I starts with ! or ?, means boolean argument
            % '!name' equivalent to 'name', false
            % '?name' equivalent to 'name', true
            % So the next element in varargin is directly a new field name !
            
            if (regexp(param, '^\?'))
                
                param = param(2:end);
                val = true;
                
            elseif (regexp(param, '^\!'))
                
                param = param(2:end);
                val = false;
                
                
            elseif (regexp(param, '\s|='))
                
                i = regexp(param, '\s|=', 'once');
                val = param(i+1:end);
                param = param(1:i-1);
                
            else
                
                if (k == numel(varargin))
                    
                    error('Parameter nb %d (%s) is not followed by a value', k, param);
                    
                end
                
                k = k + 1;
                val = varargin{k};
                
            end
            
            % Split field nameand put results in a cell array of strings
            
            Param = regexp(param, '\.', 'split');
            
            % In older versions of matlab, regexp(..., 'split')
            % returns the empty list if the seperator is absent
            % Here is the fix for that
            
            if (isempty(Param))
                
                Param = {param};
                
            end
            
            % Update Args by adding the pair (param, value) into it
            
            % Create temparry cell aray O the contains the existing tree structure
            % of Arg that converges to the parameter to be added
            
            O = cell(1 + numel(Param), 1);
            
            d = 1;
            O{d} = Args;
            
            while (d <= numel(Param) && isfield(O{d}, Param{d}))
                
                
                O{d+1} = O{d}.(Param{d});
                
                
                d = d + 1;
                
            end
            
            if (d > numel(Param))
                
                % Parameter appears more than once
                warning ('readopt', ...
                    'The argument ''%s'' appears more that once (the last passed value is retained!)', param);
                
            end
            
            % Complete tree structure from top
            
            D = numel(Param);
            O{D+1} = val;
            
            while (D > d)
                
                try
                    
                    O{D} = struct();
                    O{D}.(Param{D}) = O{D+1};
                    
                catch
                    
                    % Param{D} is not a valid field name ...
                    
                    error('Unable to set field ''%s'', is this a valid field name : ''%s'' ??', param, Param{D});
                    
                end
                
                D = D - 1;
                
            end
            
            % Merge added and existing tree structures
            
            
            if (~isstruct(O{D}))
                
                O{D} = struct();
                O{D}.(Param{D}) = O{D+1};
                D = D-1;
                
            end
            
            while (D > 0)
                
                O{D}.(Param{D}) = O{D+1};
                D = D - 1;
                
            end
            
            % Write O{1} into Args
            
            Args = O{1};
            
            
        end
        
        k = k + 1;
        
    end
end

% Split Args

[tArgs, ntArgs] = fill_args(args, Args);

end

function [Args, AddArgs] = fill_args(T, args)

if (~isstruct(T))
    
    AddArgs = [];
    
    if (iscell(T) && ~iscell(args))
        
        Args = {args};
        
    else
        
        Args = args;
        
    end
    
    return;
    
end

if (~isstruct(args))
    
    Args = T;
    AddArgs = args ;
    
    return;
end

f = fieldnames(T);
g = fieldnames(args);
h = setdiff(g, f);

Args = struct();
AddArgs = struct();

for k = 1:numel(f)
    
    if (isfield(args, f{k}))
        
        [Args.(f{k}), addArgs] = fill_args(T.(f{k}), args.(f{k}));
        
        if (~isempty(addArgs))
            
            AddArgs.(f{k}) = addArgs;
            
        end
        
    else
        
        Args.(f{k}) = T.(f{k});
        
    end
    
end

for k = 1:numel(h)
    
    AddArgs.(h{k}) = args.(h{k});
    
end
end

function param = findmatch(param, argf)

% Find The first alphabetical character in string param
i = regexp(param, '[a-zA-Z]', 'once');

% Find the first spacing character or equal sign in input string param
j = regexp(param, '\s|=', 'once');
if (isempty(j))
    j = length(param);
else
    j = j - 1;
end

if (isempty(i) || j < i)
    
    error ('Bad shortcut parameter ''%s'', not of the form ''-(!|?)<id>'' !!', param);
    
end

% Find an entry in argf that starts with param(i:j)

mparam = argf(strmatch(param(i:j), argf));

if (numel(mparam) > 1)
    error('Ambiguous ''shortcut'' parameter ''%s'', %d matches found !', param, numel(mparam));
elseif (isempty(mparam))
    error('No match found for ''shortcut'' parameter ''%s'' !', param);
end

% Transform param into the fully determined parameter name

param = [param(2:i-1), mparam{1}, param(j+1:end)];

end
