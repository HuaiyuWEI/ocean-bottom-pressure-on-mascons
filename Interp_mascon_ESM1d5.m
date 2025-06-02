% Huaiyu Wei
% This script is to interpolate ACCESS-ESM1.5 model outputs from its native 
% grid to the mascon grid used by JPL.
% 

clear all
close all
addpath(genpath('D:\OneDrive - University of California\MATLAB toolboxs'));
addpath(genpath('D:\OneDrive - University of California\MATLAB Codes\MOC\CMIP6_2025'));


%% load basin mask

BasinMasksFN = 'D:\OneDrive - University of California\MATLAB Codes\MOC\CMIP6\BasinMask\BasinMasks_CM4.mat';
load(BasinMasksFN)

pcolor(X,Y,MaskAtl)
colorbar
shading flat


%% load the mascon grid

grdfile = 'D:\OneDrive - University of California\PythonWork\mascon\GRCTellus.JPL.200204_202503.GLO.RL06.3M.MSCNv04CRI.nc';
nc=netcdf(grdfile);
ncdisp(grdfile)

lon_mascon = nc{'lon'}(:);
lat_mascon = nc{'lat'}(:);
[lon_mascon,lat_mascon] = meshgrid(lon_mascon,lat_mascon);
mascon_ID = nc{'mascon_ID'}(:);
land_mask = nc{'land_mask'}(:);

% remove mascons in polar region, and on land
mascon_ID(or(lat_mascon<-75,lat_mascon>64.5)) = nan;
mascon_ID_uniq = unique(mascon_ID(~isnan(mascon_ID)));
for i = 1:length(mascon_ID_uniq)
    ind =  (mascon_ID == mascon_ID_uniq(i));
    if(mean(land_mask(ind)) ~=0)
    mascon_ID_uniq(i) = nan;
    end

end


% label all mascon
mascon_ID_uniq = mascon_ID_uniq(~isnan(mascon_ID_uniq));
Nmascon = length(mascon_ID_uniq);

% find the lat lon for each labeled mascon
lon_mascon_bound1 = zeros(Nmascon,1);
lon_mascon_bound2 = zeros(Nmascon,1);
lat_mascon_bound1 = zeros(Nmascon,1);
lat_mascon_bound2 = zeros(Nmascon,1);
for i = 1:Nmascon
    ind =  (mascon_ID == mascon_ID_uniq(i));
    lon_mascon_bound1(i) = min(lon_mascon(ind)) - 0.25;
    lon_mascon_bound2(i) = max(lon_mascon(ind)) + 0.25;
    lat_mascon_bound1(i) = min(lat_mascon(ind)) - 0.25;
    lat_mascon_bound2(i) = max(lat_mascon(ind)) + 0.25;
end

lon_mascon_center = 0.5*(lon_mascon_bound1+lon_mascon_bound2);
lat_mascon_center = 0.5*(lat_mascon_bound1+lat_mascon_bound2);



% identify the basin of each mascon
Atlantic_id = nan(Nmascon, 1);  % Preallocate output
lon_mascon_center_180 = lon_mascon_center;
lon_mascon_center(lon_mascon_center>180) = lon_mascon_center(lon_mascon_center>180)-360;
for i = 1:Nmascon
    dist2 = (X - lon_mascon_center(i)).^2 + (Y - lat_mascon_center(i)).^2;
    % Find the index of the closest point
    [~, ind] = min(dist2(:));
    [row, col] = ind2sub(size(MaskAtl), ind);
    % Assign the basin mask value
    Atlantic_id(i) = MaskAtl(row, col);
end
Atlantic_id( Atlantic_id==2)=nan;
scatter(lon_mascon_center, lat_mascon_center, 15,  Atlantic_id, 'filled');
xlabel('Longitude')
ylabel('Latitude')
colorbar;
%%
ind_str = 1;
ind_end =5;
% DataPath_ESM1d5 = 'E:\Data_CMIP6\ACCESS_historical\';
% DataPath_ESM1d5 = 'E:\Data_CMIP6\ACCESS_SSP585\';
% DataPath_ESM1d5 = 'E:\Data_CMIP6\ACCESS_SSP245\';
% DataPath_ESM1d5 = 'E:\Data_CMIP6\ACCESS_SSP370\';
DataPath_ESM1d5 = 'E:\Data_CMIP6\ACCESS_SSP126\';
% DataPath_ESM1d5 = 'E:\Data_CMIP6\ACCESS_SSP126\2100-2300\';
% DataPath_ESM1d5 = 'E:\Data_CMIP6\ACCESS_SSP585\2100-2300\';
%% load input variables from ACCESS-ESM1.5
InputNames = {'pbo'};
NInputs = length(InputNames);
cd(DataPath_ESM1d5)

% Get list of 'pbo' files to determine ensemble members
files = dir('pbo_Omon_ACCESS-ESM1-5*.nc');
file_names = {files.name};
% Extract number after "_r" (e.g., r2, r11, etc.)
r_numbers = zeros(size(file_names));
for i = 1:length(file_names)
    % Match '_r<num>i' using regular expression
    tokens = regexp(file_names{i}, '_r(\d+)i', 'tokens');
    if ~isempty(tokens)
        r_numbers(i) = str2double(tokens{1}{1});
    end
end
% Sort the files based on the r_numbers
[~, sort_idx] = sort(r_numbers);
filelist = file_names(sort_idx);
% Display sorted list
% disp(filelist')



for i = ind_str:ind_end%length(filelist)
    file_base = filelist{i};
    fprintf('Processing realization file: %s\n', file_base);

    % Extract realization number
    tokens = regexp(file_base, 'r(\d+)i1p1f1', 'tokens');
    realization = '000';  % default if no match
    if ~isempty(tokens)
        realization = tokens{1}{1};
    end

    Input_vars = [];
    %%% Load input variables and apply missing value mask
    for var_ind = 1:NInputs
        InputName = InputNames{var_ind};
        file_name = fullfile(DataPath_ESM1d5, [InputName, file_base(4:end)]);
        nc=netcdf(file_name);
        Var=nc{InputName}(:);
        fullvalue = ncreadatt(file_name,InputName,'_FillValue');
        Var(Var==fullvalue)=nan;
        Input_vars(var_ind,:,:,:) =Var;
    end

    % Reorder to [lon x lat x time x var]
    Input_vars = permute(Input_vars, [4, 3, 2, 1]);
    Input_vars = [Input_vars(281:360,:,:,:); Input_vars(1:280,:,:,:)];

    %%% Load coordinates and time (only once)
    if i == ind_str
        sample_file = fullfile(DataPath_ESM1d5, ['pbo', file_base(4:end)]);
        Input_time = ncread(sample_file, 'time');
        Nsamps_ESM1d5 = length(Input_time);
        Lat_ESM = ncread(sample_file, 'latitude');
        Lon_ESM = ncread(sample_file, 'longitude');
        Lon_ESM = [Lon_ESM(281:360,:); Lon_ESM(1:280,:)];
        Lat_ESM = [Lat_ESM(281:360,:); Lat_ESM(1:280,:)];
    end


   
  % Average OBP from ACCESS's grid onto the mascons 
  Input_vars_mascon = zeros(Nmascon,Nsamps_ESM1d5);
   for i = 1:Nmascon
    ind =  Lon_ESM>=lon_mascon_bound1(i) & Lon_ESM<lon_mascon_bound2(i) ...
         & Lat_ESM>=lat_mascon_bound1(i) & Lat_ESM<lat_mascon_bound2(i);
     % Find the row and column subscripts of these indices
    [row, col] = find(ind);
    Input_vars_mascon(i, :) = squeeze(nanmean(Input_vars(row, col,:),[1 2]));
    end

 % obtain OBP_mascon(lat,lon)
Input_vars_mascon_lonlat = nan*zeros(size(mascon_ID,1),size(mascon_ID,2),Nsamps_ESM1d5);
for i = 1:Nmascon
ind =  (mascon_ID == mascon_ID_uniq(i));
[row, col] = find(ind);
Input_vars_mascon_lonlat(row, col, :) = repmat(reshape(Input_vars_mascon(i, :),[1,1,Nsamps_ESM1d5]),[length(row),length(col),1]);
end
  
% have a look on the data
figure
pcolor(lon_mascon,lat_mascon,Input_vars_mascon_lonlat(:,:,88))
shading flat
colorbar

    %%% Save result
    save(['Mascon_OBP_r' realization '.mat'], ...
        'lon_mascon','lat_mascon', 'Input_vars_mascon_lonlat','Atlantic_id', ...
        'Input_vars_mascon', 'lon_mascon_center','lat_mascon_center', ...
        'lon_mascon_bound1','lon_mascon_bound2','lat_mascon_bound1','lat_mascon_bound2','-v7.3');


    %%% Clean up for next loop
    clear Input_vars Input_vars_gr F
end


