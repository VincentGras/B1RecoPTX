function [ima, nlin, ncol] = mosaic(ima, nlin, ncol)

% 3d-image -> 2d-mosaic

if (nargin <= 2)
    ncol = 0;
end

if (nargin <= 1)
    nlin = 0;
end

if (isscalar(ima)) 
    n = ima;
else
    n = numel(ima)/size(ima, 1)/size(ima,2);
end

if (nlin)
    ncol = ceil(n/nlin);
elseif (ncol)
    nlin = ceil(n/ncol);
else
    nlin = 1:sqrt(n);
    nlin = nlin(find(mod(n./nlin,1)==0, 1, 'last'));
    ncol = ceil(n/nlin);
    if (ncol /nlin > 3)    
        nlin = ceil(sqrt(n));
        ncol = ceil(n/nlin);
    end
end

if (isscalar(ima))
    return;
end

if (nlin * ncol > n)
    ima = padarray(ima, [0,0,nlin*ncol-n], NaN, 'post');
end

%ima = permute(abs(ima), [2 1 3]);
ima = reshape(ima, [size(ima,1), size(ima, 2), ncol, nlin]);
ima = permute(ima, [1 4 2 3]);
ima = reshape(ima, size(ima, 1) * nlin, size(ima, 3) * ncol);
