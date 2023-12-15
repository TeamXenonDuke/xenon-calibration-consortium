%Script to plot views from Siemens data like
% (meant to replicate good old pfile_read_sort from the GE era)

% Version history:
%{ 
 -2/14/17 extended to allow for multi-slice data
 -7/7/17 expand header fields read and deal with multi-channel stuff
 -1/6/17 update to deal with UVA Headers
 -3/23/21 attempt to correct ref_voltage and num_coils
 -3/30/21 corrected ref_voltage output for more sequences
 -4/6/21 attempt to correct FFT scale factor output for some sequences
 -6/15/21 update to read Siemens Vida (software XA20A) files
 -6/28/21 option to read bonus calibration FIDs on end of Dixon seq
 -7/1/21 auto identifies if bonus spectra on dixon, prints key bonus cali
   vars
 -8/24/21 standardizes field name "sWiPMemBlock" (capital P) -- requires
   RenameField function
 -9/23/21 fix bug when reading UVA dixon '2007' files
 -10/26/21 BD display peak and RMS signals of first frame for QA, remove
   Lorentz barrier fit option, remove ref amp calculation for bonus cali
   spectra
- 12/19/21 BD Add display of pulse length and offset frequency

%}

% Required functions:
%{
-mapVBVD.m
-RenameField.m
-ReadBonusSpectra.m
%}

% Classes:
%{
-twix_map_obj.m
%}

% Bastiaan Driehuys, Aryil Bechtel
%%
clc; clear all; close all;
plot_phases=1;                 % plot real imaginary and magnitude data (useful for spectra)
first_frames=5;                % number of initial fids to plot typically 5-20
ymax = 1;                      % if not set to one, then set yscale of first figure to this value
duke_spectro = 0 ;             % set to 1 only for basic duke spectroscopy (Soher) sequence to calculate n_fids properly. 
data_fraction = 1;             % plot smaller fractions of data to keep giant data sets from crashing the computer

%% first read in the file and create the twix object
[file, path] = uigetfile('*.*', 'Select file');  % starts in current dir 
file_with_path = strcat(path, file);  % join path and filename to open
twix_obj = mapVBVD(file_with_path); 
filename = file(15:end-4); % return filename back to something short
filename = strrep(filename, '_', '-'); % replace underscores to avoid subscript problems

%if twix_obj is 1xn (n>1) cell, choose which cell
if size(twix_obj,2)>1
    optAns = input('noise or image? N/I:','s');
    if strcmp(optAns,'N')
        twix_obj = twix_obj{1};
        twix_obj = RenameField(twix_obj,'noise','image');
    elseif strcmp(optAns,'I')
        twix_obj = twix_obj{2};
    end
end
%% if twix obj contains field sWipMemBlock, change to sWiPMemBlock (capital P)

if isfield(twix_obj.hdr.MeasYaps,'sWipMemBlock')
    temp = RenameField(twix_obj.hdr.MeasYaps,'sWipMemBlock','sWiPMemBlock');
    twix_obj.hdr.MeasYaps = temp;
end

%% if Dixon, determine if there exist bonus spectra

%if adFree{6} is integer and adFree{11} contains info, assume bonus spectra
if and(contains(filename,'Dixon'),~contains(filename,'BHUTE')) 
    if and(isfield(twix_obj.hdr.MeasYaps.sWiPMemBlock,'adFree'),...
            length(twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree)>5)
        numDisSpect = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{6};
        numGasSpect = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{11}; 
        if and(mod(numDisSpect,1)==0,~isempty(numGasSpect))
            containsBonus = 1;
            numSpect = numDisSpect + numGasSpect;
        else
            containsBonus = 0;
        end
    else 
        containsBonus = 0;
    end
else
    containsBonus = 0;
end

%% dealing with fields that may or may not be in header for all sequences
if(~isfield(twix_obj.hdr.Meas, 'flTransRefAmpl'))
      header.imagetwix_obj.hdr.Meas.flTransRefAmpl = 0;
end

if(~isfield(twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}, 'flAmplitude'))
      twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude = 0;
end

if(~isfield(twix_obj.hdr.Meas, 'GradDelayTime'))
      twix_obj.hdr.Meas.GradDelayTime = 0;
end

if(~isfield(twix_obj.hdr.Config, 'NCha'))
      twix_obj.hdr.Config.NCha = 0;
end

if(~isfield(twix_obj.hdr.Config, 'nCha'))  % if nCha also isn't a field, then just one channel
      twix_obj.hdr.Config.nCha = 1;
end


if(~isfield(twix_obj.hdr.Config, 'NSlc'))
      twix_obj.hdr.Config.NSlc = 0;
end

if(~isfield(twix_obj.hdr.Config, 'NoOfFourierPartitions'))  % if reading 3D data
      twix_obj.hdr.Config.NoOfFourierPartitions = 1;
end

if(~isfield(twix_obj.hdr.Meas, 'lGain'))
      twix_obj.hdr.Meas.lGain = 0;
end

%% read key header variables
Seq_name = twix_obj.hdr.Config.SequenceDescription;
weight = twix_obj.hdr.Dicom.flUsedPatientWeight; 
freq = twix_obj.hdr.Dicom.lFrequency;  % seems to work in all sequences
%pulse_length = twix_obj.hdr.Phoenix.sWipMemBlock.adFree{1}; % in us
freq_offset = twix_obj.hdr.Phoenix.sWipMemBlock.alFree{5} % dissolved offset frequency

if(~isfield(twix_obj.hdr.Meas, 'alTE(1)')) % doesn't seem to exist in UVA data
    TE = twix_obj.hdr.Phoenix.alTE{1}; % works for UVA Dixon
    %TE=450; % just hardwire it
else
    TE = twix_obj.hdr.Meas.alTE(1);
end %if    
TR = twix_obj.hdr.Config.TR;
dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};  % appears most robust and generic

% get or calculate reference voltage from twix object
if and(contains(filename,'Dixon'),~contains(filename,'BHUTE')) 
    if contains(filename,'BW799')
        if size(twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE,2)>1
            rf_amp2 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,2}.flAmplitude; %assuming 2nd rfamp in sTXSPEC is the correct one
            ref_voltage = 600*rf_amp2/193.2;
        else
            ref_voltage = 0;
        end
    elseif contains(filename,'2008')
        if size(twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree,2)>2
            ref_voltage = twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{1,3};
        else
            ref_voltage = 0;
        end
    else
        if isfield(twix_obj.hdr.MeasYaps.sTXSPEC,'aRFPULSE')
            rf_amp1 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude;
            ref_voltage = 600*rf_amp1/193.2; % estimate reference voltage assuming 0.69ms sinc pulse
        else
            ref_voltage = 0;
        end
    end
elseif contains(filename,'2DGRE')
    if size(twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree,2)>3
        ref_voltage = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{1,4};
    else
        ref_voltage = 0;
    end
elseif and(contains(filename,'ventilation'),str2double(filename(end))<=9 && str2double(filename(end))>=0)
    if size(twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree,2)>2
        ref_voltage = twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{1,3};
    else
        ref_voltage = 0;
    end
elseif contains(filename,'ventilation')
    if size(twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE,2)>1
        rf_amp2 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,2}.flAmplitude; %assuming 2nd rfamp in sTXSPEC is the correct one
        ref_voltage = 600*rf_amp2/193.2;
    else
        ref_voltage = 0;
    end
elseif or(contains( filename, 'cali' ),contains(filename,'Cali'))
    if isfield(twix_obj.hdr.MeasYaps.sTXSPEC,'aRFPULSE')
        rf_amp1 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude;
        ref_voltage = 600*rf_amp1/193.2;
    else
        ref_voltage = 0;
    end
else
    if isfield(twix_obj.hdr.Spice,'TransmitterReferenceAmplitude')
        ref_voltage = twix_obj.hdr.Spice.TransmitterReferenceAmplitude;
    else
        ref_voltage = 0;
    end
end

% get fft scale factor
if(isfield(twix_obj.hdr.MeasYaps,'asCoilSelectMeas')) % not in UVA data
    fft_factor = twix_obj.hdr.MeasYaps.asCoilSelectMeas{1,1}.aFFT_SCALE{1,1}.flFactor;
elseif(isfield(twix_obj.hdr.MeasYaps,'sCoilSelectMeas'))
    fft_factor = twix_obj.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.aFFT_SCALE{1,1}.flFactor;
else
    fft_factor = 0;
end %if

rec_gain = twix_obj.hdr.Meas.lGain; % receiver gain, high if 1, low if 0
grad_delay= twix_obj.hdr.Meas.GradDelayTime;  % a factor in radial images that seems to vary
num_coils = twix_obj.image.NCha; % number of coil channels
num_slices = twix_obj.hdr.Config.NSlc;  % number of slices
num_partitions = twix_obj.hdr.Config.NoOfFourierPartitions; % allow for 3D data 
if num_partitions > num_slices  % it's 3D  
    num_slices = num_partitions; % just call 'em slices and move on
end %if
num_pts = twix_obj.image.sqzSize(1);  % num_pts is always first in array. Second can be coils or FIDs
num_fids = twix_obj.hdr.Config.NoOfFourierLines; % works for imaging sequences, but not basic Duke spectroscopy
if duke_spectro == 1
    num_fids = twix_obj.hdr.Config.NRepMeas; % this works for basic duke spectroscopy sequence. The above contains number of points for spectro
end % if

% for multislice GRE can use Config.ImageLines, NImagLins, PhaseEncodingLines, RawLn, NoOfFourierLines 
% for radial can use Config.NLinMeas, NProj, NoOfFourierLines, 
% doesn't yet deal with averaging

%% print out key header variable values
if containsBonus==1
    fprintf('\n %s\n','Key variables for Dixon seq (excluding bonus seqs)');
end
fprintf('\n %s\n',file_with_path);  % show the full path for patient ID info
fprintf('\n');
fprintf('\tFile = %s\n',filename');
fprintf('\tSequence = %s\n',Seq_name');
fprintf('\tWeight = %0.1f kg \n',weight);
fprintf('\tTE=%0.0f us\n',TE);
fprintf('\tTR=%0.0f us\n',TR(1));
fprintf('\tFrequency=%0.0f Hz\n',freq);
fprintf('\treference voltage=%0.1f V\n',ref_voltage);
%fprintf('\tpulse_length = %0.0f us\n',pulse_length);
fprintf('\tfreq_offset = %0.0f Hz\n',freq_offset);
fprintf('\tfft_scale_factor=%0.3f\n',fft_factor);
fprintf('\tdwell_time=%0.0f ns\n',dwell_time);
fprintf('\tnum_coils=%0.0f\n',num_coils);
fprintf('\tnum_slices=%0.0f\n',num_slices);
fprintf('\tnPts=%0.0f\n',num_pts);
fprintf('\tnFrames=%0.0f\n',num_fids);

%% create the actual FID data and deal with dimensions and slices
twix_obj.image.flagIgnoreSeg = true; %Add this flag to ignore the extra radial dimensions. Doesn't hurt other files.
theFID = squeeze(double(twix_obj.image.unsorted()));  % use for giant arrays, just lays slices in there.

if ndims(theFID) == 3 % this happens for multi-slice or multi-channel data
    if num_coils > 1 % multi coil data presents as num_pts x num_chans x (num_slicesxnum_fids?
        iCoil = 1; % just look at one coil at a time for now. Change this to look at different coil
        theFID2 = squeeze(theFID(:,1,:)); % eliminates coil dimension and makes matrix of nptsx(nslicesxnframes) compatible with rest of data
        fprintf('\nPlotting only coil %g\n',iCoil);
    else              % multi-slice data - when does this happen?
        num_fids1 = num_fids*num_slices;
        theFID2=reshape(theFID,num_pts, num_fids1); % just line up all the fids from multi slice
    end % if
    
else  % just simple 1-D data
    theFID2=theFID; % just keep theFID as the name
end %if

if data_fraction ~= 1
    stop_point = round(data_fraction*size(theFID,2));
    theFID=theFID(:,1:stop_point);
    fprintf('\nPlotting only first %2.1f of data\n',data_fraction);
end %if 

if containsBonus==1 %extract bonus spectra if they exist
    theFIDlength = max(size(theFID));
    theFID(:,theFIDlength-numSpect:end) = []; %remove poor res bonus FIDs
    theFID2(:,theFIDlength-numSpect:end) = [];
    [bonusCali] = readBonusSpectra(twix_obj); %get full res bonus FIDs
end

%% if bonus cali FIDs, get key variables and print

if containsBonus==1
    % get spec ratios and plot spectral fits
    nDis = numDisSpect; % # dissolved FIDs
    nSkip = 1; % # dissolved FIDs to skip
    nAvg = numDisSpect - 2; % # dissolved FIDs to avg
    bonusCali = double(bonusCali);
    nPts = size(bonusCali,1);
    bonusDwell = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{10}*0.5*1e-6;
    calRefVolt = twix_obj.hdr.MeasYaps.sWiPMemBlock.alFree{3};
    disTR = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{7} * 1e-3;
    gasTR = twix_obj.hdr.MeasYaps.sWiPMemBlock.adFree{12} * 1e-3;
    
    fprintf('\n %s\n','Key variables for bonus cali seq');
    fprintf('\tdisTR=%0.0f us\n',disTR*1e6);
    fprintf('\tgasTR=%0.0f us\n',gasTR*1e6);
    fprintf('\tdwell_time=%0.0f ns\n',bonusDwell*1e9);
    fprintf('\tnPts=%0.0f\n',nPts);
    fprintf('\tnumDisFIDs=%0.0f\n',numDisSpect);
    fprintf('\tnumGasFIDs=%0.0f\n',numGasSpect);
end

%% separate magnitude plots of various sections of data 
h1=figure;   % first number of fids requested
if first_frames>num_fids
    first_frames=num_fids;
    sprintf('\nNumber of plot frames exceeded. Setting to %g\n',first_frames);
end 
plot(abs(theFID2(1:first_frames*num_pts)));    % plot only first_frames
axis tight;
xlabel('points')
ylabel('magnitude')
a=sprintf('%s First %0.0f fids',filename,first_frames);
if ymax ~= 1
    ylim([0 ymax])
end 
title(a)
box on
fprintf('Peak signal in first frames = %3.3e\n',max(abs(theFID2(1:first_frames*num_pts))));
fprintf('RMS signal of first frames = %3.3e\n',rms(abs(theFID2(1:first_frames*num_pts))));


h2=figure;   % plot all fids sequentially
if containsBonus==1 %append bonus FIDs if they exist
    allFIDs = [abs(theFID2(1:end)) abs(bonusCali(1:end))];
else
    allFIDs = abs(theFID2(1:end));
end
plot(allFIDs); 
hold on;
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All FIDs Sequential',filename);
title(a)
box on

h3=figure;   % plot all Dixon fids overlapping to look for noise
plot(abs(theFID2)); 
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All Dixon FIDS Overlapping',filename);
if ymax ~= 1
    ylim([0 ymax])
end %if
title(a)
box on

%if bonus FIDs, plot dis and gas phase cali FIDs overlapping
if containsBonus==1
    %plot dis phase cali FIDs overlapping
    h5 = figure;
    bonusDis = bonusCali(:,1:numDisSpect);
    plot(abs(bonusDis));
    xlabel('points')
    ylabel('magnitude')
    a=sprintf('%s All Bonus Dissolved FIDS Overlapping',filename);
    title(a)
    box on;

    %plot gas phase cali FIDs overlapping
    h6 = figure;
    bonusGas = bonusCali(:,(numDisSpect+1):end);
    plot(abs(bonusGas));
    xlabel('points')
    ylabel('magnitude')
    a=sprintf('%s All Bonus Gas FIDS Overlapping',filename);
    title(a)
    box on;
end

% plot to get both real and imaginary channels
if plot_phases==1
    h4=figure;
    ax1=subplot(311);
    plot(real(theFID2(1:end)));
    grid on;
    xlabel('points');
    ylabel('Real Signal');
    a=sprintf('%s - Phase Sensitive Display',filename);
    title(a)
    
    ax2=subplot(312);
    plot(imag(theFID2(1:end)));
    grid on;
    xlabel('points');
    ylabel('Imaginary Signal');
    
    ax3=subplot(313);
    %plot(abs(theFID_array(1:end)));
    plot(abs(theFID2(1:end)));
    grid on;
    xlabel('points');
    ylabel('Magnitude');
    
    linkaxes([ax1,ax2,ax3],'x'); % just link x axes since magnitude is different

end %if

figure(h1); % bring figure 1 back to the front of the pack
