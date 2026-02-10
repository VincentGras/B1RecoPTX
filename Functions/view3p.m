function view3p(map, varargin)

if (nnz(imag(map(:)))>0)
    warning('Taking magnitude');
    map = abs(map);
end

DimX = size(map, 1);
DimY = size(map, 2);
DimZ = size(map, 3);


% read options    

opt.layout = {1,3};
opt.cscale = [];
opt.mosaic = {};
opt.s = [0,0,0];
opt.o = 'c';
opt.ctitle = ''; 
opt.axis = 'image';
opt.newfig = [];
opt.subplot = 1;

opt = readopt(opt, varargin{:});

if (strcmpi(opt.o, 'c') || strcmpi(opt.o, 'center') || strcmpi(opt.o, 'centre'))
    opt.s = round([size(map,1), size(map,2), size(map,3)]/2 + opt.s);
end

opt.s = max(opt.s,1);
opt.s = min(opt.s, [size(map,1), size(map,2), size(map,3)]);

if (isempty(opt.cscale))
    opt.cscale = double([min(map(:)), max(map(:))]);
    if (opt.cscale(2) <= opt.cscale(1))
        opt.cscale(2) = opt.cscale(1)+1;
    end
end



mapSagittal = flip(permute(squeeze(map(opt.s(1),:,:,:)), [2 1 3]), 1);
mapCoronal = flip(permute(squeeze(map(:,opt.s(2),:,:)), [2 1 3]), 1);
mapAxial = permute(squeeze(map(:,:,opt.s(3),:)), [2 1 3]);

if (isempty(opt.newfig))
    opt.newfig = opt.subplot == 1;
end

if (opt.newfig)
    figure();
end

subplot(opt.layout{:},opt.subplot);
data = mosaic(mapSagittal, opt.mosaic{:});
imagesc(data, 'AlphaData', ~isnan(data)); 
axis (opt.axis); 
caxis(opt.cscale);
title('sagittal'); 
xlabel('A \rightarrow P'); 
ylabel('F \leftarrow H'); 
set(gca, 'XTick', []);
set(gca, 'YTick', []);
colorbar;
subplot(opt.layout{:},opt.subplot+1);
data = mosaic(mapCoronal, opt.mosaic{:}); 
imagesc(data, 'AlphaData', ~isnan(data)); 
axis (opt.axis); 
caxis(opt.cscale);
title('coronal');
ylabel('F \leftarrow H');
xlabel('R \rightarrow L'); 
set(gca, 'XTick', []);
set(gca, 'YTick', []);
colorbar;
subplot(opt.layout{:},opt.subplot+2);
data = mosaic(mapAxial, opt.mosaic{:});
imagesc(data, 'AlphaData', ~isnan(data)); 
axis (opt.axis); 
caxis(opt.cscale); 
colormap(jet);
title('axial');
xlabel('R \rightarrow L');
ylabel('P \leftarrow A'); 
set(gca, 'XTick', []);
set(gca, 'YTick', []);
colorbar;
