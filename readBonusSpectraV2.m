function [raw,numGasSpect,numDisSpect] = readBonusSpectraV2(source)
% readBonusSpectraSimple - works with Twix or ISMRMRD
% 
% Usage:
%   raw = readBonusSpectraSimple(twix_obj);               % Twix
%   raw = readBonusSpectraSimple('file.h5',10,20);        % MRD file (10 dis + 20 gas)
%
% Returns:
%   raw : [FIDlength x numSpect] complex
[~, base, ext] = fileparts(source);

if isstruct(source) || strcmpi(ext, '.dat')
 twix_obj = mapVBVD(source); 
if isfield(twix_obj.hdr.MeasYaps,'sWipMemBlock')
    spectReso = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{9}; % read in spect sample points
    numDisSpect = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{6}; % find # of bonus dis spectra
    numGasSpect = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{11}; % find # of bonus gas spectra 
elseif isfield(twix_obj.hdr.MeasYaps,'sWiPMemBlock')
    spectReso = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{9};
    numDisSpect = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{6};
    numGasSpect = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{11}; 
end
numSpect = numDisSpect + numGasSpect;
obj = twix_obj.image;
imgReso = obj.dataSize(1)/2; % read in image sample points

% extract parameters for spectrum extraction
mem = obj.memPos; %start byte # for each FID (this is an array) 
szScanHeader = obj.freadInfo.szScanHeader; %# bytes in header to skip
readSize     = obj.freadInfo.sz;
readSize(2) = readSize(2)-imgReso*2+spectReso*2; % # elements (not bytes) in file designated for each bonus FID
readCut      = obj.freadInfo.cut; %readCut(1) gives # designated elements not part of each bonus FID 

% open .dat file
fid = obj.fileopen();

%initialize output matrix
rawLength = readSize(2) - readCut(1) + 1;
raw = zeros(rawLength,numSpect);

% skip scan header and extract data
for i=1:numSpect
    fseek(fid,mem(i+numel(mem)-numSpect) + szScanHeader,'bof'); %set file reading positionz
    raw1 = fread(fid, readSize, 'float=>single').';
    raw1 = complex(raw1(:,1), raw1(:,2));
    raw(:,i) = raw1(readCut(1):end,:); %chop off readCut(1) elements not part of each bonus FID
end


elseif ischar(source) || any(strcmpi(ext, {'.h5','.mrd'}));

    import ismrmrd.*

    dset = ismrmrd.Dataset(char(source),'dataset');
    nAcq = dset.getNumberOfAcquisitions();

fids   = {};
labels = [];

for i = 1:nAcq
    acq = dset.readAcquisition(i);

    % Keep only bonus spectra
    if acq.head.measurement_uid == 1
        if iscell(acq.data)
            chData = cell2mat(acq.data);   % [nChannels x nSamples]
        else
            chData = acq.data;
        end
        fids{end+1} = chData; %#ok<AGROW>
        labels(end+1) = double(acq.head.idx.contrast); %#ok<AGROW>
    end
end

% Build matrix
nSamp = max(cellfun(@numel, fids));
numSpect = numel(fids);
raw = zeros(nSamp, numSpect);
for k = 1:numSpect
    fk = fids{k};
    raw(1:numel(fk), k) = fk;
end

% Counts
numGasSpect = sum(labels == 1);
numDisSpect = sum(labels == 2);



else
    error('Unsupported input type');
end
end
