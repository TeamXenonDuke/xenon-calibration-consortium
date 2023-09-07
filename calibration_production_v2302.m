% Script to plot and fit calibration FIDs

% Required functions:
%{
-mapVBVD.m
-RenameField.m
-NMR_TimeFit.m
-NMR_TimeFit_v.m
-readRawCali.m
%}

% Classes:
%{
-twix_map_obj.m
%}


clc;clear all;close all;

%% Define number of calibration, spectral fit and display toggles,number to skip, number of gas, target flip

% user control of  outputs
disp_fids = 1; % set to 1 to display first (noise) frame and first FID
save_csv =1; % export derived variables as csv file

% user control of skipping and averaging
seconds2skip = 2; %number of seconds to skip in breath hold
seconds2avg = 1;
nDis = 500; % Assume consortium standard
FlipTarget = 20; % target flip angle from calibration

% frequency guesses in ppm for dissolved phase fits
rbc_freq_ppm = 217.2;
mem_freq_ppm = 197.7;
gas_freq_ppm = 0;

%% define tolerances for warnings
freq_tol = 100; % frequency adjustment tolerance
deltaPhase1_tol = 90; % if calibration phase>90 use minTE
TrueRefScale_tol = 1.33;
SNR_tol = 50; % Ideal is 100 for RBC SNR, but allow 50

%% first read in the file and create the twix object
[file, path] = uigetfile('*.*', 'Select file'); % starts in current dir
file_with_path = strcat(path, file); % join path and filename to open
[~,filename,ext] = fileparts(file_with_path); % get file extension
cali_struct = readRawCali(file_with_path);
if strcmp(ext,'.dat') 
   filename = file(15:end-4); % return filename back to something short
else
   filename = [filename,'_cali'];
end   
cal_vars_export = [filename,'.csv']; % generate the filename for data comparison
filename = strrep(filename, '_', '-'); % replace underscores to avoid subscript problems
file_loc = regexp(path, filesep, 'split');
file_loc = file_loc{end-1};

%% Prepare variables
% extract favorite variables from the data struct
weight = cali_struct.weight;
te = cali_struct.te;
tr = cali_struct.tr;
dwell_time = cali_struct.dwell_time;
freq = cali_struct.freq;
xeFreqMHz = cali_struct.xeFreqMHz;
theFID = cali_struct.data;
nFids = size(theFID, 2);
nCal = nFids-nDis; % assume remaining FIDS past dissolved are cal
nPts = size(theFID, 1);
VRef = cali_struct.vref;
scanDateStr = cali_struct.scan_date;
rf_excitation_ppm = cali_struct.rf_excitation_ppm;

% print out key variable values
fprintf('\n\n');
fprintf('\tFile Location = %s\n', file_loc);
fprintf('\tName = %s\n', cali_struct.seq_name);
fprintf('\tScan Date = %s\n', scanDateStr');
fprintf('\tWeight = %0.1f kg \n', weight);
fprintf('\tTE = %0.0f us\n', te);
fprintf('\tTR = %0.0f us\n', tr(1)); % some sequences now have an array of TRs for bonus spectra
fprintf('\tFrequency = %0.0f Hz\n', freq);
fprintf('\tRF Excitation Offset = %0.0f ppm\n', rf_excitation_ppm);
fprintf('\tReference Voltage = %0.1f V\n', VRef);
fprintf('\tDwell Time = %0.0f ns\n', dwell_time*1e9);
fprintf('\tnPts = %0.0f\n', nPts);
fprintf('\tnFrames = %0.0f\n', nFids);
fprintf('\tnDissolved Frames = %0.0f\n', nDis);
fprintf('\tnGas Frames = %0.0f\n\n', nCal);

%% parse out the various fids - dissolved, gas, and calibration
calData = theFID(:, end-nCal+1:end); % the data left for flip angle calculations
gasData = theFID(:, end-nCal+1); % the first gas frame for frequency calculations

tr_s = tr(1) * 1e-6; %tr in seconds
t_tr = tr_s * (1:nFids);
nSkip = round(seconds2skip/tr_s);
nAvg = round(seconds2avg/tr_s);
t = dwell_time * (0:(nPts - 1))';
disData = theFID(:, nSkip:nSkip+nAvg); % just the dissolved data of interest
disData_avg = mean(disData,2); % average the dissolved data

%% Fit gas Spectrum
fprintf('Analysis of Gas FID\n');
gasfitObj = NMR_TimeFit(gasData, t, 1e-4, -84, 30, 0, 0, 10000);
gasfitObj.fitTimeDomainSignal();
h1 = figure('Name', 'Gas Phase Analysis');
gasfitObj.plotTimeAndSpectralFit;
xlim([0, 0.01]);
gasfitObj.describe();

%% calculate target transmit receive frequency from gas fit
freq_target = freq + gasfitObj.freq(1);
fprintf('\nFrequency_target = %8.0f Hz\n', freq_target); % report new target frequency
if abs(freq_target-freq) > freq_tol
    fprintf(2, 'WARNING! Frequency adjust exceeds tolerances; Check system\n');
end

%% Perform Flip Angle Calibration on Gas Phase Signals
flipCalAmps = max(abs(calData));

% calculate flip angle
fitfunct = @(coefs, xdata)coefs(1) * cos(coefs(2)).^(xdata - 1) + coefs(3); % cos theta decay
guess(1) = max(flipCalAmps);
guess(2) = 20 * pi / 180; % just guess 10 degrees
guess(3) = 0; % with absolute values there will be a baseline offset per Matt Wilmering 4/12/21

xdata = 1:length(flipCalAmps);
ydata = flipCalAmps;

fitoptions = optimoptions('lsqcurvefit', 'Display', 'off');
[fitparams, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(fitfunct, guess, xdata, ydata, [], [], fitoptions);
ci = nlparci(fitparams, residual, jacobian); % returns 95% conf intervals on fitparams by default
param_err = fitparams - ci(:, 1)';
flip_angle = abs(fitparams(2)*180/pi);
flip_err = param_err(2) * 180 / pi;

h2=figure('Name', 'Flip Angle Calibration'); % plot up the calibration data
plot(xdata, ydata, 'bo', 'MarkerFaceColor', 'b');
hold on;
plot(xdata, fitfunct(fitparams, xdata), '-r');
legend('Acquired', 'Fit');
xlabel('Frame Number');
ylabel('Magnitude');
a = sprintf('Flip Cal: V_{ref} = %0.1f yields flip = %0.1f�%0.1f%c', VRef, flip_angle, flip_err, char(176));
title(a); % note degree symbol is char(176)
pause(1); % had to insert pause to avoid occasional plotting conflicts - don't know why

%% Provide New Reference Amplitude and Warnings
VRefScaleFactor = FlipTarget / flip_angle; % How much should Ref Voltage Be Scaled?
if isnan(VRef)
   hint_refV = '\bf Reference Voltage: \rmNot found';
   fprintf(2,'WARNING! Reference Voltage not found.\n');
else   
   fprintf(['For ', num2str(VRef), 'V calibration, True_Ref = %3.0f V\n'], VRef*VRefScaleFactor); 
   fprintf('(Estimated Weight Reference Voltage = %3.0f V)\n', 395.3+1.1431*weight); % Ari Bechtel estimate, new fit
   hint_refV = ['\bf Reference Voltage: \rm', num2str(VRef*VRefScaleFactor, 3), 'V'];
end   
if VRefScaleFactor > TrueRefScale_tol
    fprintf(2, 'WARNING! Excessive voltage calibration scale factor; Check system\n');
end 

%% Fit dissolved Spectrum - Make this last so user sees it
fprintf('\nAnalysis of Averaged Dissolved FID, skipping %0.0f\n', nSkip);
fprintf('Using Voigt Lineshape to fit Membrane\n');

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

% first widths lorenzian, 2nd are gaussian
disfitObj = NMR_TimeFit_v(disData_avg, t, area_guess, freq_guess, fwhmL_guess, fwhmG_guess, phase_guess, [],[]);
disfitObj = disfitObj.fitTimeDomainSignal();
h3 = figure('Name', 'Dissolved Phase Analysis'); 
disfitObj.plotTimeAndSpectralFit;
% update axis limits to see gas signal
ax = findall(h3, 'type', 'axes'); % find all the axes in figure (ax1 defined in NMR_TimeFit)
set(ax(1),'xlim',[0 .01]); % narrow gas plot limits for meaningful time
set(ax(5),'xlim',[-8000 2000])  % expand plot limits to show gas peak

disp('            Area     Freq (Hz)   Linewidths(Hz)   Phase(degrees)');
peakName = {'     RBC:', 'Membrane:', '     Gas:'};
for iComp = 1:length(disfitObj.area)
    disp([peakName{iComp}, '  ', ...
        sprintf('%8.3e', disfitObj.area(iComp)), ' ', ...
        sprintf('%+8.2f', disfitObj.freq(iComp)), '  ', ...
        sprintf('%8.2f', abs(disfitObj.fwhm(iComp))), '  ', ...
        sprintf('%8.2f', abs(disfitObj.fwhmG(iComp))), '  ', ...
        sprintf('%+9.2f', disfitObj.phase(iComp))]);
end

%% Calculate and report TE90 and Warnings
deltaPhase = disfitObj.phase(2) - disfitObj.phase(1); % RBC-Membrane phase diff in cal spectrum
deltaPhase = mod(abs(deltaPhase), 180); % deal with wrap around, but also negative phase
deltaF = abs(disfitObj.freq(2)-disfitObj.freq(1)); % absolute RBC-membrane freq difference
deltaTe90 = (90 - deltaPhase) / (360 * deltaF); % how far off are we?
te90 = (te + deltaTe90 * 1e6)/1000; % in usec

% Hint the actual TE90 for usage
if te90 < 0.445 
   hint_te90 = '\color[rgb]{1,0,0} WARNING! Low TE90; Use 0.45ms.';
elseif  te90 >= 0.505
   hint_te90 = '\color[rgb]{1,0,0} WARNING! High TE90; Use 0.50ms.';
else
   hint_te90 = ['\color[rgb]{0,0,0} Use ',num2str(te90, 2), 'ms.'];
end    
fprintf('\nTE90 = %3.2f ms\n', te90); % report TE90 in ms
fprintf(2,[hint_te90(20:end),' \n']);

%% Calculate and report dissolved phase spectral SNR
pause(0.1);
disData_tail = disData(nPts/2+1:end,:); % just take the second half of FIDs
mean_disData_tail = mean(disData_tail,2);
diff = mean_disData_tail - disData_tail; % array of difference fids
diff = diff - mean (diff); % subtract off any DC bias
noiseFrame = diff(:,1); % take first frame as representative noise frame
noiseDis = std(real(diff));
SNR_frames_d = disfitObj.area'./noiseDis; % SNR of each frame
SNRsnf_d = mean(SNR_frames_d,2)*sqrt(nAvg); % simulated noise frame SNR
noise_mean = mean(noiseDis,2); % save mean noise frame for gas calculation

fprintf('\nSNR Check for dissolved and gas peaks...\n'); 
fprintf('SNR of dissolved spectra (Dai method):\n');
fprintf('SNR for RBC peak = %3.1f \n', SNRsnf_d(1));
fprintf('SNR for Membrane peak = %3.1f \n', SNRsnf_d(2));
fprintf('SNR for Gas peak = %3.1f \n', SNRsnf_d(3));

% Calculate and report SNR of analyzed gas FID using the Dai method
SNR_gas_frame = gasfitObj.area'./noise_mean; % use noise calculated from dissolved
fprintf('SNR for dedicated gas peak = %3.1f\n',SNR_gas_frame); 

%% Quantify RBC:Membrane ratio
RbcMemRatio = disfitObj.area(1) / disfitObj.area(2);
fprintf('\nRbcMemRatio = %3.3f\n', RbcMemRatio); %

%% Quantify ammount of off resonance excitation
GasDisRatio = disfitObj.area(3) / sum(disfitObj.area(1:2));
fprintf('GasDisRatio = %3.3f \n', GasDisRatio);

%% Save some parameters for comparing if requested
if save_csv == 1
    folderPath = path; % save the csv file to the same directory as the dat file
    fileName = cal_vars_export; % name of the file
    outputFilePath = fullfile(folderPath, fileName); % full path

    % define a struct to save bonusSpetraCalibration calculations
    cal_cali.FreqTarget = freq_target;
    cal_cali.GasSNR = SNR_gas_frame;
    cal_cali.TrueRefV = VRef*VRefScaleFactor; %reference voltage
    cal_cali.TE90 = te90;
    cal_cali.RbcMemRatio = RbcMemRatio;
    cal_cali.GasDisRatio = GasDisRatio;
    cal_cali.RBCShift = (disfitObj.freq(1) - disfitObj.freq(3)) / freq_target * 1e6; %frequency difference between RBC and gas in ppm
    cal_cali.RBCSNR = SNRsnf_d(1);
    cal_cali.MemSNR = SNRsnf_d(2);

    %write and save calculated variables to file
    dataTable_cali = struct2table(cal_cali);
    writetable(dataTable_cali, outputFilePath, 'WriteVariableNames', true);
end 

%% Display Representative Noise Frame and First FID
% if this is done before all the fitting it messes up axis scaling
if disp_fids == 1 % display representative fids only if desired
    h4 = figure('Name', 'Noise Frame Analysis'); % plot first frame for time-domain noise check
    plot(abs(theFID(1:nPts)));  % plot only the first_frame
    hold on;
    axis tight;
    ymax=max(abs(theFID(1:nPts)));
    if ymax<5e-5 % Siemens scale
        ylim([0 1e-5]) % standard Duke noise scale
    end 
    plot(nPts/2+1:nPts,abs(noiseFrame)); % add the simulated noise (2nd half) frame
    xlabel('points')
    ylabel('magnitude')
    box on
    title('Noise Frame Analysis')
    legend('Dedicated','Simulated')
    rmsNoise=[rms(abs(theFID(1:nPts))) rms(abs(noiseFrame))];
    fprintf('\nRMS noise signal (dedicated) = %3.3e \n',rmsNoise(1));
    fprintf('RMS noise signal (simulated) = %3.3e \n',rmsNoise(2));
    h5 = figure('Name', 'First FID'); % plot up first representative FID for time-domain noise check
    plot(real(theFID(:,2)), '-b');
    hold on;
    plot(imag(theFID(:,2)), '-r');
    legend('Real', 'Imag');
    xlabel('Points');
    ylabel('Real Signal');
    a = sprintf('First Signal Frame');
    title(a);
end

set(0, 'currentfigure', h3); shg; % bring the dissolved figure to front and show it

%% PJN - Creat Pop-up box of parameters that need to be set in real time
CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';
textstr = {['\bf Gas Center Frequency: \rm', num2str(freq_target, 8), ' Hz']; ' ';...
    hint_refV;' ';...
    ['\bf TE_{90}: \rm', num2str(te90, 2), ' ms']; ' ';...
    hint_te90;};
outbox = msgbox(textstr, 'Imaging Parameters', CreateStruct);

th = findall(outbox, 'Type', 'Text'); %get handle to text within msgbox
th.FontSize = 18; %Change the font size

deltaWidth = sum(th.Extent([1, 3])) - outbox.Position(3) + th.Extent(1);
deltaHeight = sum(th.Extent([2, 4])) - outbox.Position(4) + 10;
outbox.Position([3, 4]) = outbox.Position([3, 4]) + [deltaWidth, deltaHeight];

outbox.Resize = 'on';

% Move the message box away from the plot
figPosition = outbox.Position; % Get the position of the current figure
msgboxX = figPosition(1) + figPosition(3); % Set the X position for the message box
msgboxY = figPosition(2); % Set the Y position for the message box
outbox.Position(1:2) = [msgboxX, msgboxY]; % Update the position of the message box
