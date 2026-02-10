function out = vect2map(msk, dat, fillval)

% Convert vectorized image into 3D
% see also: mat2vect

if (nargin < 3)
    fillval = 0;
end

if (iscell(msk))
    msk = cat(4, msk{:});
end

msk = logical(msk);

msk_sz = size(msk);

if (ismatrix(msk))
    msk_sz(3) = 1;
end

%out = cell(size(varargin));
%for i = 1:numel(out)
%out{i} = zeros(size(msk));
%dat = squeeze(varargin{i});

dat_sz = size(dat);

if (ismatrix(dat))
    if (size(dat, 1) ~= nnz(msk(:)))
        dat = dat.';
        dat_sz = size(dat);
    end    
else
    dat = dat(:,:);
end

assert(nnz(msk) == size(dat, 1), 'data size is not compatible with mask size');

nv = size(dat, 2);

if (nv == 1)
    
    out = zeros(size(msk), class(dat));
    out(:) = fillval;
    out(msk) = dat;
    
else
    
    out = cell(nv, 1);
    
    for j = 1:nv
        out{j} = zeros(size(msk), class(dat));
        out{j}(:) = fillval;
        out{j}(msk) = dat(:, j);
    end
    out = cat(min(ndims(msk),3)+1, out{:});
    
end

out = reshape(out, [msk_sz, dat_sz(2:end)]);

%end

%out = cat(4, out{:});

