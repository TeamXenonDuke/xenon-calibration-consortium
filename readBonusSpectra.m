function [raw] = readBonusSpectra(twix_obj)

%{
modification of read_twix_bonus_spect_ZW.m (Ziyi 17.08.14)
modified by Aryil Bechtel 2021

intended for use with new Dixon sequence:
    -10 appended dissolved FIDs
    -20 appended dedicated gas FIDs

input: radial Dixon twix_obj with the numSpect extra FIDs
output: (FIDlength x numSpect matrix) of FIDs
%}

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
    fseek(fid,mem(i+numel(mem)-numSpect) + szScanHeader,'bof'); %set file reading position
    raw1 = fread(fid, readSize, 'float=>single').';
    raw1 = complex(raw1(:,1), raw1(:,2));
    raw(:,i) = raw1(readCut(1):end,:); %chop off readCut(1) elements not part of each bonus FID
end

end

