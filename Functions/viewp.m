function viewp(p, map, varargin)

% usage:
% ima = view3dp(p, map, varargin)
% where options are :
% cscale 1x2 double
% title  string
% mosaic cell containing mosaic parameters (see mosaic.m)
% s      1x3 double position in volume

if (nnz(imag(map(:)))>0)
    warning('Taking magnitude');
    map = abs(map);
end

opt.cscale = [min(map(:)), max(map(:))];
if (opt.cscale(2)<=opt.cscale(1))
    opt.cscale(1) = opt.cscale(1) - 0.5;
    opt.cscale(2) = opt.cscale(2) + 0.5;
end
opt.title = '';
opt.mosaic = {};
opt.s = 0;
opt.o = 'c';
opt.ctitle = '';
opt.axis = 'image';
opt.colorbar = true;
%opt.display = true;
opt.layout = {}; %{1,1};
opt.newfig = [];
opt.subplot = 1;

opt = readopt(opt, varargin{:});

if (strcmpi(opt.o, 'c') || strcmpi(opt.o, 'center') || strcmpi(opt.o, 'centre'))
    opt.s = round(size(map, p)/2 + opt.s);
end

opt.s = max(opt.s,1);
opt.s = min(opt.s, size(map, p));

if (ischar(opt.s) && strcmp(opt.s, '*'))
    opt.s = 1:size(map, p);
end

if (ndims(map)>4)
    warning ('map dimensionnality is > 4');
    map = map(:,:,:,:);
end

DimX = size(map, 1);
DimY = size(map, 2);
DimZ = size(map, 3);

if (p == 1)
    ima = flip(permute(map(opt.s,:,:,:), [3 2 4 1]), 1);
    xl = 'A \rightarrow P';
    yl = 'F \leftarrow H';
elseif (p == 2)
    ima = flip(permute(map(:,opt.s,:,:), [3 1 4 2]), 1);
    yl = 'F \leftarrow H';
    xl = 'R \rightarrow L';
    
elseif (p == 3)
    ima = permute(map(:,:,opt.s,:), [2 1 4 3]);
    xl = 'R \rightarrow L';
    yl = 'P \leftarrow A';
    
end

data = mosaic(ima(:, :, :), opt.mosaic{:});

% if (~opt.display)
%     return;
% end

if (isempty(opt.newfig))
    opt.newfig = opt.subplot == 1;
end

if (opt.newfig)
    figure();
end

if (~isempty(opt.layout))
    subplot(opt.layout{:},opt.subplot);
end
imagesc(data, 'AlphaData', ~isnan(data));
axis (opt.axis);
caxis(opt.cscale);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
xlabel(xl);
ylabel(yl);
if (opt.colorbar)
    hcb = colorbar;
    ylabel(hcb, opt.ctitle);
end
title (opt.title);
