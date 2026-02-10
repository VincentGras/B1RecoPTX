%% Test resulting FA distributions for the encoding scheme

addpath .\Functions\

db = '.\Avanti2_B1maps_Agar.map.mat'; % file with simulated B1+ fields (Nx x Ny x Nz x Nch, T/V)
gyr = 42577000;
dataB1 = load(db);
data.M = dataB1.Mask; % mask
data.b1 = map2vect(dataB1.Mask,  dataB1.B1map).' * (2*pi) * gyr; % (transform to vector form Nch x Nvox, rad/s/V)

paramfile='.\Schemes\Sphere_vcc03_Acond2_Bcond2.xflparam.txt'; % file with encoding scheme


fid = fopen(paramfile,'r');
numLines = 4;
for ii = 1:numLines
    line = fgetl(fid); 
    if (~isempty(regexpi(line, 'Vref')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        Vref=str2num(a);
    end

    if (~isempty(regexpi(line, 'alpha_nominal')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        alpha_nom=str2num(a);
    end

    if (~isempty(regexpi(line, 'beta_nominal')))
        a=regexp(line, '(\d*)', 'match', 'once' );
        beta_nom=str2num(a);
    end
end
fclose(fid);


[alphaMat, betaMat] = readXFLParametersFromFile(paramfile);

RFPulseIntegralAlpha = alpha_nom/180 * Vref * 1e-3 *180/pi; % Volt.s
alphaMat = RFPulseIntegralAlpha * alphaMat; 

RFPulseIntegralBeta = beta_nom/180 * Vref * 1e-3 *180/pi; % Volt.s
betaMat = RFPulseIntegralBeta * betaMat;

sat = alphaMat(:,any(alphaMat,1));

FA_alpha=abs(sat'*data.b1);
FA_beta=abs(betaMat'*data.b1);


mean_alpha=mean(abs(FA_alpha(:)));
mean_alpha_ind=mean(abs(FA_alpha),2).';
max_alpha=max(FA_alpha(:));
min_alpha=min(FA_alpha(:));
min_max_alpha=min(max(FA_alpha,[],1));

fprintf( strcat("mean alpha for each scheme = ", repmat('%.3f ', 1, numel(mean_alpha_ind)),"\n"), mean_alpha_ind );
fprintf("min alpha =  %.2f \n", min_alpha );
fprintf("max alpha =  %.2f \n", max_alpha );
fprintf("min of max alpha over the LCCs =  %.2f \n", min_max_alpha );

view3p(applyMask(data.M, vect2map(data.M,  FA_alpha)),'cscale', [0 180]);
annotation('textbox',[0.4 0.7 0.1 0.1],'String', '\alpha','FontSize',20,'EdgeColor','none') %04 038 039

view3p(applyMask(data.M, vect2map(data.M,  max(FA_alpha,[],1))),'cscale', [0 180]); 
annotation('textbox',[0.4 0.7 0.1 0.1],'String', 'max \alpha','FontSize',20,'EdgeColor','none') %04 038 039

figure, 
for i=1:size(FA_alpha,1)
    subplot(2,4,i)
    histogram(FA_alpha(i,:),[0:2:200]); title('\alpha distribution'); xlabel('FA, °'); ylabel('Nb of voxels')
end

Max_FA= max(FA_alpha,[],1);
figure, histogram(abs(Max_FA(:))); title('max \alpha distribution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



mean_beta_ind=mean(abs(FA_beta),2).';
max_beta=max(FA_beta(:));
min_beta=min(FA_beta(:));
min_max_beta=min(max(FA_beta,[],1));

fprintf( strcat("mean beta for each scheme = ", repmat('%.3f ', 1, numel(mean_beta_ind)),"\n"), mean_beta_ind );
fprintf("min beta =  %.2f \n", min_beta );
fprintf("max beta =  %.2f \n", max_beta );
fprintf("min of max beta over the LCCs =  %.2f \n", min_max_beta );

view3p(applyMask(data.M, vect2map(data.M,  FA_beta)),'cscale', [0 20]);
annotation('textbox',[0.4 0.7 0.1 0.1],'String', '\beta','FontSize',20,'EdgeColor','none') %04 038 039

view3p(applyMask(data.M, vect2map(data.M,  max(FA_beta,[],1))),'cscale', [0 20]); 
annotation('textbox',[0.4 0.7 0.1 0.1],'String', 'max \beta','FontSize',20,'EdgeColor','none') %04 038 039

figure,
for i=1:size(FA_beta,1)
    subplot(9,2,i);
    histogram(abs(FA_beta(i,:)),[0:0.1:15]); title('\beta distribution'); xlabel('FA, °'); ylabel('Nb of voxels')
end

Max_Beta= max(FA_beta,[],1);
figure, histogram(FA_beta(:)); title('max \beta distribution')

