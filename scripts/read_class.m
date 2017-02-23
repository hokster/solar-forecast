function data = read_class(folder,type,varargin)
% folder is the path to the .nc file
% type is the data to be read out: 
%   'data' is the cloud radiance data
%   'lat' is the latitude data
%   'lon' is the longitude data
%   'time' is the time data for the image.
% varargin is an optional vector of dates in the format used in the .nc
% files

% read in the files --- all if varargin is empty, and just the entries in
% varargin if nonempty.
if isempty(varargin)
    files = dir(fullfile(folder,'goes13*.nc'));    
else
    for i=1:length(varargin)
        files(i) = dir(fullfile(folder,sprintf('goes13*%s*.nc',varargin{i})));            
    end
end

% read in data from the collected files
for i = 1:length(files)
    if strcmp(type,'data')
        data(:,:,i) = ncread(fullfile(folder,files(i).name),'data');
    elseif strcmp(type,'lat')
        data(:,:,i) = ncread(fullfile(folder,files(i).name),'lat');
    elseif strcmp(type,'lon')
        data(:,:,i) = ncread(fullfile(folder,files(i).name),'lon');
    elseif strcmp(type,'time')
        data(i) = ncread(fullfile(folder,files(i).name),'time');
    end
end