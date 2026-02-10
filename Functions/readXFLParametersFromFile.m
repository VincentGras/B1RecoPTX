function [rfcoeffsat, rfcoeffexc, ref] = readXFLParametersFromFile(fid)

if (ischar(fid))
    
    fname = fid;
    fid = fopen(fname);
    assert(fid>0, 'cannot read file ''%s''', fname);
    
    try
        
        [rfcoeffsat, rfcoeffexc, ref] = readXFLParametersFromFile(fid);
        fclose(fid);
    catch err
        fclose(fid);
        rethrow(err);
    end
    return;
    
end

assert(gotoXFLParameters(fid), 'cannot find XFL parameters');

rfcoeffsat = cell(0,1);
rfcoeffexc = cell(0,1);
satpulse = zeros(1,0);

cnt = 0;

while (true)
    
    l = nextline(fid);
    
    if (~ischar(l))
        break;
    end
    
    if (~isempty(regexpi(l, '^SATURATE_RFSAT')))
       
        cnt = cnt + 1;
        c = readcoeff(l);
        rfcoeffsat{cnt} = c;
        l = nextline(fid);
        assert(ischar(l), 'file incompleted');
        assert(~isempty(regexpi(l, '^SATURATE_RFEXC')), 'after a SATURATE_RFSAT line, expecting SATURATE_RFEXC line');
        c = readcoeff(l);
        rfcoeffexc{cnt} = c;
        satpulse(cnt) = 1; % scan #cnt is a "sat pulse"
        
    elseif (~isempty(regexpi(l, '^RELATIVE_RFEXC')))
        
        cnt = cnt + 1;
        c = readcoeff(l);
        rfcoeffexc{cnt} = c;
        rfcoeffsat{cnt} = 0*c;
        satpulse(cnt) = 0; % scan #cnt is not a "sat pulse"
    else
        
        error('cannot interprete line ''%s''', l);
        
    end    
    
end

% transform into matrices #lines = #TX channels & #col = #of scans
rfcoeffsat = cat(2, rfcoeffsat{:});
rfcoeffexc = cat(2, rfcoeffexc{:});

% define ref as follows :
% ref(i) = -1 if scan #i is not as "sat" scan
% ref(i) = 0 if scan #i is a "sat" scan , but sat pulse has zero amplitude (all(rfcoeffsat(:,i) == 0)) 
% ref(i) = j>0 if scan #i is a "sat" scan with a sat pulse with non-zero
%              amplitude for at least one TX channel. In this case j return
%              the scan index for the corresponding reference: 
%                 (sat(j)==true && all(rfcoeff(:,i) == rfcoeff(:,j)))

ref = zeros(1, size(rfcoeffsat, 2));
voidsat = all(rfcoeffsat==0,1);

for i = 1:numel(ref)
    
    if (satpulse(i))
        
        if (any(rfcoeffsat(:,i)))
            
            d = max(abs(bsxfun(@minus, rfcoeffexc(:,i), rfcoeffexc)), [], 1);
            j = find(d < 1e-9 & satpulse & voidsat);
            if (numel(j) == 0)
                warning ('cannot find reference scan for satscan #%d', i);
                j = 0;
            elseif (numel(j)>1)
                warning ('multiple reference scans for satscan #%d (%s)', i, num2str(j));
                if (ismember(i-1, j))
                    j = i-1; % prefer the scan right before
                else
                    j = j(1);
                end
            end
            ref(i) = j;
        else 
            
            ref(i) = 0;
        
        end
        
    else
        
            ref(i) = -1;
        
    end
end

function l = nextline(fid, skipcomment)

skipline = @(l) (l(1) == '#');

if (nargin >= 2 && ~skipcomment)
    skipline = @(l) false;
end

while (true)
    
    l = fgetl(fid);
    
    if (~ischar(l))        
        return;
    end
    
    l = strip(l);
    
    if (isempty(l))
        continue;
    end
    
    if (skipline(l))
        continue;
    end
    
    return;
    
end

function c = readcoeff(l)

[~, j] = regexpi(l, '^[a-z_]+\s*=?');
c = sscanf(l(j+1:end), '%f', inf);
c = c(1:2:end) .* exp(1i * pi/180 * c(2:2:end));

function success = gotoXFLParameters(fid)

success = false;

while (true)
    
    l = nextline(fid, false);
    
    if (~ischar(l))
        break;
    end
    
    if (~isempty(regexpi(l, 'XFLparameters', 'once')))
        
        success = true;
        break;
        
    end
end
    
    



