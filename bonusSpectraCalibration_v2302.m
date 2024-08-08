%{
adapted by Aryil Bechtel (2021) from demo_calibration_duke_UVA_B1.m
updated by Yi Zheng (2023)

script to run calibration analysis on bonus calibration FIDs.
currently intended for use with 30 bonus calibration FIDs.

dissolved FID spectral fit:
    -skips 1st FID
    -averages (# diss FIDs)-2
gas FID spectral fit:
    -uses 2nd FID
%}

clc; clear all; close all;

%% Value setttings for spectral fit and display toggles, tolerances for warnings, target flip
Voigt = 1;% 0 = Membrane Lorentzian, 1 = Membrane Voigt fitting
          % Requires NMR_fit_v, NMR_mix_v, and NMR_TimeFit_v for Voigt=1
nSkip = 1; % # dissolved FIDs to skip - bd note: why?

% frequency guesses in ppm for dissolved phase fits
rbc_freq_ppm = 217.2;
mem_freq_ppm = 197.7;
gas_freq_ppm = 0;

FlipTarget = 20; % target flip angle from calibration
gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla

%% Read in the file and create the twix object
[file, path] = uigetfile('*.*', 'Select file');  % starts in current dir 
file_with_path = strcat(path, file);  % join path and filename to open
twix_obj = mapVBVD(file_with_path); 
filename = file(15:end-4); % return filename back to something short
cal_vars_export = [filename,'.csv']; % generate the filename for data comparison
filename = strrep(filename, '_', '-'); % replace underscores to avoid subscript problems
file_loc = regexp(path,filesep,'split');
file_loc = file_loc{end-1};

%% Check if Dixon and determine if bonus spectra exist

%if adFree{6} is integer and adFree{11} contains info, assume bonus spectra
if contains(filename, 'Dixon', 'IgnoreCase', true) && ~contains(filename, 'BHUTE') && ~contains(filename, 'proton', 'IgnoreCase', true)
    if isfield(twix_obj.hdr.MeasYaps,'sWiPMemBlock')
       twix_obj.hdr.MeasYaps.sWipMemBlock = twix_obj.hdr.MeasYaps.sWiPMemBlock; %replace the old name
       twix_obj.hdr.MeasYaps = rmfield(twix_obj.hdr.MeasYaps,'sWiPMemBlock');
    end    
    numDisSpect = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{6};
    numGasSpect = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{11}; 

    if and(mod(numDisSpect,1)==0,~isempty(numGasSpect))
        containsBonus = 1;
        numSpect = numDisSpect + numGasSpect;
    else
        error('ERROR: INCOMPATIBLE INPUT FILE, must be a radial Dixon + cali sequence');
    end
else
    error('ERROR: INCOMPATIBLE INPUT FILE, must be a radial Dixon + cali sequence');
end

%% Extract bonus FIDs and key variables
[bonusCali] = readBonusSpectra(twix_obj); %get full res bonus FIDs
Seq_name = twix_obj.hdr.Config.SequenceDescription;
weight = twix_obj.hdr.Dicom.flUsedPatientWeight;
freq = twix_obj.hdr.Dicom.lFrequency; 
te = twix_obj.hdr.Phoenix.alTE{1}; %for now, assuming this is correct TE for bonus FIDs
dixDwell = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}*1e-9; %dixon dwell time in seconds

% get spec ratios and plot spectral fits
nDis = numDisSpect; % # dissolved FIDs
nAvg = numDisSpect - 2; % # dissolved FIDs to avg

bonusCali = double(bonusCali);
nPts = size(bonusCali,1); % # pts in each bonus FID
nFids = numSpect;
nCal = numGasSpect;

bonusDwell = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{10}*0.5*1e-6; %dwell time in seconds, divide 2 bc oversampling
calRefVolt = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{3};
disTR = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{7} * 1e-3; %dissolved TR in seconds
gasTR = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{12} * 1e-3; %gas TR in seconds 

% get scan date
scanDate = twix_obj.hdr.Phoenix.tReferenceImage0; 
scanDate = strsplit(scanDate,'.');
scanDate = scanDate{end};
scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];

% print out key variable values
fprintf('\n\n');
fprintf('\tName = %s\n',Seq_name');
fprintf('\tScan Date = %s\n',scanDateStr');
fprintf('\tWeight = %0.1f kg \n',weight);
fprintf('\tTE = %0.0f us\n',te);
fprintf('\tdisTR = %0.0f us\n',disTR*1e6);
fprintf('\tgasTR = %0.0f us\n',gasTR*1e6);
fprintf('\tFrequency = %0.0f Hz\n',freq);
fprintf('\tReference Voltage = %0.1f V\n',calRefVolt);
fprintf('\tDwell Time = %0.0f ns\n',bonusDwell*1e9);
fprintf('\tnPts = %0.0f\n',nPts);
fprintf('\tnumDisFIDs=%0.0f\n',numDisSpect);
fprintf('\tnumGasFIDs=%0.0f\n',numGasSpect);

% parse out the various fids - dissolved, gas, and calibration
disData = bonusCali(:,1:nDis); % all dissolved data
gasData = bonusCali(:,end-nCal+2);
if (nSkip+1+nAvg) > nDis
    fprintf('\n Requested averages exceeds available; using %0.0f\n',nDis);
    disData1 = bonusCali(:,nSkip+1:nDis); % all available dissolved data with skipping
else
    disData1 = bonusCali(:,nSkip+1:(nSkip+1+nAvg)); % requested averages of dissolved
end 
disData1_avg = mean(disData1,2); % average of dissolved phase
calData = bonusCali(:,end-nCal+1:end); % the data left for flip angle calculations
t = double((0:(length(disData)-1))*bonusDwell');

%% Read RF excitation frequency
% Magnetic Field Strength
mag_fstrength = twix_obj.hdr.Dicom.flMagneticFieldStrength;

% Read excitation from twix header,Dixon version
excitation = twix_obj.hdr.Phoenix.sWipMemBlock.alFree{1, 5};

% RF excitation will be in ppm, likeley either 218ppm or 208 ppm at Duke
rf_excitation_ppm = round(excitation/(gyro_ratio * mag_fstrength));

%% Fit gas Spectrum
fprintf('\nAnalysis of Gas FID\n');
gasfitObj = NMR_TimeFit(gasData,t,1e-4,-84,30,0,0,10000);
gasfitObj.fitTimeDomainSignal();
figure('Name','Gas Phase Analysis')
gasfitObj.plotTimeAndSpectralFit;
xlim([0 0.01]);
gasfitObj.describe();

%% Fit dissolved Spectrum
fprintf('\nAnalysis of Averaged Dissolved FID, no skipping\n');
xeFreqMHz = freq/1e6;

% Adjust frequency guesses based on rf_excitation then to Hz
rbc_freq_adj = rbc_freq_ppm - rf_excitation_ppm;
mem_freq_adj = mem_freq_ppm - rf_excitation_ppm;
gas_freq_adj = gas_freq_ppm - rf_excitation_ppm;

freq_guess = [rbc_freq_adj, mem_freq_adj, gas_freq_adj] * xeFreqMHz; % in Hz

%set all other initial parameter guesses
area_guess =  [1, 1, 1]; % no benefit seen in tayloring these guesses
fwhmL_guess = [8.8, 5.0, 2] * xeFreqMHz;
fwhmG_guess = [0, 6.1, 0] * xeFreqMHz;
phase_guess = [0, 0, 0]; % no benefit seen in tayloring these guesses

if Voigt == 1
   fprintf('Using Voigt Lineshape to fit Membrane\n');
   disfitObj = NMR_TimeFit_v(disData1_avg, t, area_guess, freq_guess, fwhmL_guess, fwhmG_guess, phase_guess, [],[]);% first widths lorenzian, 2nd are gauss
else
   fprintf('Using Lorentzian Lineshape to fit Membrane\n'); % will work for 218 ppm excitation only
   disfitObj = NMR_TimeFit(disData1_avg,t,[1 1 1],[0 -700  -7400],[240 240 40],[0 0 0],0,length(t));
end 
disfitObj = disfitObj.fitTimeDomainSignal();
figure('Name','Dissolved Phase Analysis')
disfitObj.plotTimeAndSpectralFit;
% update axis limits to see gas signal
ax = findall(gcf, 'type', 'axes'); % find all the axes in figure
set(ax(1),'xlim',[0 .01]); % narrow plot limits for meaningful time
set(ax(5),'xlim',[-8000 2000])  % expand plot limits to show gas (ax1 defined in NMR_TimeFit)
 
disp('            Area     Freq (Hz)   Linewidths(Hz)   Phase(degrees)');
peakName = {'     RBC:','Membrane:','     Gas:'};
if Voigt==1 %print Gaussian fwhm
    for iComp = 1:length(disfitObj.area)
        disp([peakName{iComp} '  '...
            sprintf('%8.3e',disfitObj.area(iComp)) ' ' ...
            sprintf('%+8.2f',disfitObj.freq(iComp))  '  ' ...
            sprintf('%8.2f',abs(disfitObj.fwhm(iComp))) '  ' ...
            sprintf('%8.2f',abs(disfitObj.fwhmG(iComp))) '  ' ...
            sprintf('%+9.2f',disfitObj.phase(iComp))]);
    end
else %no Gaussian fwhm to print
    for iComp = 1:length(disfitObj.area)
        disp([peakName{iComp} '  '...
            sprintf('%8.3e',disfitObj.area(iComp)) ' ' ...
            sprintf('%+8.2f',disfitObj.freq(iComp))  '  ' ...
            sprintf('%8.2f',abs(disfitObj.fwhm(iComp))) '  ' ...
            sprintf('%+9.2f',disfitObj.phase(iComp))]);
    end
end 


%% Calculate target transmit receive frequency from gas fit
freq_target = freq + gasfitObj.freq(1);
fprintf('\nFrequency_target = %8.0f Hz\n',freq_target); % report new target frequency

%% Calculate and report TE90 NEW
deltaPhase = disfitObj.phase(2)-disfitObj.phase(1); % RBC-membrane phase diff in cal spectrum
deltaPhase = mod(abs(deltaPhase),180); % deal with wrap around, but also negative phase
deltaF = abs(disfitObj.freq(2)-disfitObj.freq(1)); % absolute RBC-membrane freq difference
deltaTe90= (90-deltaPhase)/(360*deltaF); % how far off are we?
te90 =te + deltaTe90*1e6; % in usec
fprintf('TE90 = %3.2f ms\n',te90/1000); % report TE90 in ms
if abs(deltaPhase)> 90
    fprintf(2,'WARNING! Phi_cal = %3.0f%c!; Use min TE!\n',deltaPhase,char(176));
end %if

%% Flip angle calibration 
% Find Amplitudes faster by using max amplitudes from each FID in calibration
flipCalAmps = max(abs(calData));

% calculate flip angle
fitfunct = @(coefs,xdata)coefs(1)*cos(coefs(2)).^(xdata-1)+coefs(3);   % cos theta decay
guess(1)=max(flipCalAmps);
guess(2)=20*pi/180;       % just guess 20 degrees
guess(3)=0;    % allow a constant baseline per Matt Wilmering 4/12/21

xdata=1:length(flipCalAmps);
ydata = flipCalAmps;

fitoptions = optimoptions('lsqcurvefit','Display','off');
[fitparams,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqcurvefit(fitfunct,guess,xdata,ydata,[],[],fitoptions);
ci = nlparci(fitparams,residual,jacobian);  % returns 95% conf intervals on fitparams by default
param_err=fitparams-ci(:,1)';
flip_angle=abs(fitparams(2)*180/pi);
flip_err=param_err(2)*180/pi;

figure('Name','Flip Angle Calibration')  % plot up the calibration data
plot(xdata,ydata,'bo','MarkerFaceColor','b');
hold on;
plot(xdata,fitfunct(fitparams,xdata),'-r');
legend('Acquired','Fit');
xlabel('Frame Number');
ylabel('Magnitude');
a=sprintf('Flip Cal: V_{ref} = %0.1f yields flip = %0.1f�%0.1f%c', calRefVolt,flip_angle,flip_err,char(176));
title(a); % note degree symbol is char(176)

%% Provide reference amplitude and warnings
calRefVoltScaleFactor=FlipTarget/flip_angle; % How much should Ref Voltage Be Scaled?
fprintf(['For ', num2str(calRefVolt), 'V calibration, True_Ref = %3.0f V\n'],calRefVolt*calRefVoltScaleFactor); % report TE90 in ms

%% Calculate and report dissolved phase spectral SNR 
% Using spectral SNR method based on Doreen Dai's ISMRM 2023 work
fid_tail = length(disData1)/2; % focus calculation 2nd half of dissolved phase
cropped_fids = disData1(fid_tail+1:end,:); 
mean_fid = mean(cropped_fids,2);
diff = mean_fid-cropped_fids; % array of difference fids
diff = diff - mean (diff); % subtract off any DC bias
noiseDis = std(real(diff));
SNR_frames_d = disfitObj.area'./noiseDis; % SNR of each frame
SNRsnf_d = mean(SNR_frames_d,2)*sqrt(nAvg); % simulated noise frame SNR
noise_mean = mean(noiseDis,2); % save mean noise frame for gas calculation

fprintf('\nSNR of peaks in dissolved spectra...\n'); % 
fprintf('SNR for RBC peak = %3.1f \n',SNRsnf_d(1)); 
fprintf('SNR for Membrane peak = %3.1f \n',SNRsnf_d(2)); 
fprintf('SNR for Gas peak = %3.1f \n',SNRsnf_d(3)); 

% Calculate and report SNR of analyzed gas FID using the Dai method
SNR_gas_frame = gasfitObj.area'./noise_mean; % use noise calculated from dissolved
fprintf('SNR for dedicated gas peak = %3.1f\n',SNR_gas_frame); 

%% Quantify RBC:Membrane ratio
RbcMemRatio = disfitObj.area(1)/disfitObj.area(2);
fprintf('\nRbcMemRatio = %3.3f\n',RbcMemRatio); % 

%% Quantify ammount of off resonance excitation
GasDisRatio = disfitObj.area(3)/sum(disfitObj.area(1:2));
fprintf('GasDisRatio = %3.3f\n',GasDisRatio); 

%% Save some parameters for comparing
folderPath = path; % save the csv file to the same directory as the dat file
fileName = cal_vars_export; % name of the file
outputFilePath = fullfile(folderPath, fileName); % full path

% define a struct to save bonusSpetraCalibration calculations 
cal_bonus.FreqTarget = freq_target;
cal_bonus.GasSNR = SNR_gas_frame;
cal_bonus.TrueRefV = calRefVolt*calRefVoltScaleFactor; %reference voltage
cal_bonus.TE90 = te90/1000;
cal_bonus.RbcMemRatio = RbcMemRatio;
cal_bonus.GasDisRatio = GasDisRatio;
cal_bonus.RBCShift = (disfitObj.freq(1) - disfitObj.freq(3)) / freq_target * 1e6; %frequency difference between RBC and gas in ppm 
cal_bonus.RBCSNR = SNRsnf_d(1);
cal_bonus.MemSNR = SNRsnf_d(2);

%write and save calculated variables to file
dataTable_bonus = struct2table(cal_bonus);
writetable(dataTable_bonus, outputFilePath, 'WriteVariableNames', true);

