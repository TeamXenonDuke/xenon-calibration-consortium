% Extract key header variables and plot raw data from MRD files.

%% Set flag arguments
clear; close all;              % clear current variables and close figures
cal_sequence=1;                % set ==1, if processing Duke calibration sequence. If !=0, assume Duke Dixon sequence
first_frames=5;                % number of initial fids to plot typically 5-20
ymax = 1;                      % if !=1, then set y-axis scale of first figure to this value
num_cal_gas_fids=15;           % number of gas FIDs at end of calibration sequence

%% Select MRD file and extract header
[file, path] = uigetfile('*.*', 'Select file'); % open UI for selecting file, starts in current dir
file_with_path = strcat(path, file); % join path and file to open

% extract MRD header
mrd_header = xml2struct(h5read(file_with_path, '/dataset/xml'));

% vendor and institution name
vendor = lower(mrd_header.ismrmrdHeader.acquisitionSystemInformation.systemVendor.Text);
institution = mrd_header.ismrmrdHeader.acquisitionSystemInformation.institutionName.Text;

%% Extract FID data

% read k-space data
all_kspace_data = h5read(file_with_path, '/dataset/data').data;

% reshape k-space data
npts = size(all_kspace_data{1}, 1) / 2;
nfids = size(all_kspace_data, 1);
fids = [];
for k = 1:nfids
    fid = all_kspace_data{k};
    i = 1;
    for j = 1:npts
        fids(j, k) = double(complex(fid(i), fid(i+1)));
        i = i + 2;
    end
end

% if data from GE scanner, take complex conjugate
if strcmp(vendor,'ge')
    fids = conj(fids);
end

% extract gas and dissolved FIDs separately
if cal_sequence==1
    % if calibration sequence, all gas FIDs are at the end
    fids_gas = fids(:, end-num_cal_gas_fids+1:end);
    fids_dis = fids(:, 1:end-num_cal_gas_fids);
else
    % if dixon sequence, gas and dissolved are interleaved
    fids_gas = fids(:, 1:2:end);
    fids_dis = fids(:, 2:2:end);
end

%% Extract key header variables

% magnetic field strength
gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
mag_strength = str2double(mrd_header.ismrmrdHeader.acquisitionSystemInformation.systemFieldStrengthu_T.Text); % field strength in T

% variables with vendor and site-specific units and locations in header (ew how disgusting)
switch vendor
    case 'ge'
        % patient ID
        patient_id = mrd_header.ismrmrdHeader.subjectInformation.patientID.Text;

        % TR
        if cal_sequence==1
            tr = str2double(mrd_header.ismrmrdHeader.sequenceParameters.TR.Text); % in s
        else
            tr = 2 * str2double(mrd_header.ismrmrdHeader.sequenceParameters.TR.Text); % in s
        end
        tr = tr * 1e3; % convert s to ms

        % TE90
        te90 = str2double(mrd_header.ismrmrdHeader.sequenceParameters.TE.Text); % in s
        te90 = te90 * 1e3; % convert s to ms

        % dwell time
        dwell_time = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1,1}.value.Text);

        % excitation frequencies
        if cal_sequence==1
            gas_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 2}.value.Text);
            dis_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
        else
            gas_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
            dis_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 4}.value.Text);
        end
        rf_excitation = (dis_freq - gas_freq)/(gyro_ratio * mag_strength);

        % institution-specific variables
        switch institution
            case "University of Iowa"
                % study date
                study_date = nan;
            case "St. Joseph's Healthcare Hamilton"
                % study date=
                study_date = mrd_header.ismrmrdHeader.studyInformation.studyDate.Text;
            case "University of Sheffield"
                % study date
                study_date = mrd_header.ismrmrdHeader.studyInformation.studyDate.Text;
        end
    case 'philips'
        % patient ID
        patient_id = nan;

        % TR
        if cal_sequence==1
            tr = str2double(mrd_header.ismrmrdHeader.sequenceParameters.TR.Text); % in ms
        else
            tr = 2 * str2double(mrd_header.ismrmrdHeader.sequenceParameters.TR.Text); % in ms
        end

        % TE90
        te90 = str2double(mrd_header.ismrmrdHeader.sequenceParameters.TE.Text); % in ms

        % dwell time
        dwell_time = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble.value.Text);

        % excitation frequencies
        gas_freq = str2double(mrd_header.ismrmrdHeader.userParameters.userParameterLong.value.Text); % in Hz
        dis_freq = nan;
        rf_excitation = str2double(mrd_header.ismrmrdHeader.userParameters.userParameterDouble.value.Text);

        % institution-specific variables
        switch institution
            case {"CCHMC", "Cincinnati"}
                % study date
                study_date = mrd_header.ismrmrdHeader.measurementInformation.seriesDate.Text;
        end
    case 'siemens'
        % patient ID
        patient_id = mrd_header.ismrmrdHeader.subjectInformation.patientID.Text;

        % TR
        tr = str2double(mrd_header.ismrmrdHeader.sequenceParameters.TR.Text); % in ms

        % TE90
        te90 = str2double(mrd_header.ismrmrdHeader.sequenceParameters.TE.Text); % in ms

        % dwell_time
        dwell_time = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1,1}.value.Text);
        if cal_sequence~=1
            dwell_time = dwell_time * 1e-1; % convert 100s of nanoseconds to us
        end

        % excitation frequencies
        if cal_sequence==1
            gas_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 2}.value.Text);
            dis_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
        else
            gas_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
            dis_freq = str2double(mrd_header.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 4}.value.Text);
        end
        rf_excitation = (dis_freq - gas_freq)/(gyro_ratio * mag_strength);

        % study date
        study_date = mrd_header.ismrmrdHeader.studyInformation.studyDate.Text;
end

% print out key header variables
fprintf('\tFile = %s\n',file');
fprintf('\tPatient ID = %s\n',patient_id');
fprintf('\tStudy date = %s\n',study_date');
fprintf('\tGas frequency = %0.0f Hz\n',gas_freq);
fprintf('\tDissolved frequency = %0.0f Hz\n',dis_freq);
fprintf('\tRF excitation = %0.0f ppm\n',rf_excitation);
fprintf('\tTR = %0.1f ms\n',tr);
fprintf('\tTE90 = %0.3f ms\n',te90);
fprintf('\tDwell time = %0.2f us\n',dwell_time);
fprintf('\tNum FIDs = %0.0f\n',nfids);
fprintf('\tNum pts in FID (npts) = %0.0f\n',npts);

%% Plot a select number of first FIDs sequentially

% initialize figure
h=figure;

% check if number of requested fids exceeds number of total fids
if first_frames>nfids
    first_frames=nfids;
    sprintf('\nNumber of plot frames exceeded. Setting to %g\n',first_frames);
end 

% plot first FIDs
plot(abs(fids(1:first_frames*npts)));
if ymax ~= 1
    ylim([0 ymax])
end 

% set title and axis labels
a=sprintf('%s First %0.0f fids',strrep(file, '_', '\_'),first_frames);
title(a)
axis tight;
xlabel('points')
ylabel('magnitude')
box on

% print max and rms of plotted FIDs
fprintf('Peak signal in first frames = %3.3e\n',max(abs(fids(1:first_frames*npts))));
fprintf('RMS signal of first frames = %3.3e\n',rms(abs(fids(1:first_frames*npts))));

%% Plot all FIDs sequentially

% initialize figure
h=figure;

% plot all FIDs sequentially
plot(abs(fids(1:end))); 

% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All FIDs Sequential', strrep(file, '_', '\_'));
title(a)
box on

%% Plot all FIDs overlapping (to look for noise)

% initialize figure
h=figure;

% plot all FIDs overlapping
plot(abs(fids)); 
if ymax ~= 1
    ylim([0 ymax])
end

% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All FIDs Overlapping', strrep(file, '_', '\_'));
title(a)
box on

%% Plot all FIDs overlapping, separating gas and dissolved

% initialize gas FID figure
h=figure;
% plot all gas FIDs overlapping
plot(abs(fids_gas)); 
if ymax ~= 1
    ylim([0 ymax])
end
% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All Gas FIDs Overlapping', strrep(file, '_', '\_'));
title(a)
box on

% initialize dissolved FID figure
h=figure;
% plot all dissolved FIDs overlapping
plot(abs(fids_dis)); 
if ymax ~= 1
    ylim([0 ymax])
end
% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All Dissolved FIDs Overlapping', strrep(file, '_', '\_'));
title(a)
box on

%% Plot real and imaginary parts of signal

% initialize figure
h=figure;
ax1=subplot(311);

% plot real part of signal
plot(real(fids(1:end)));
grid on;
xlabel('points');
ylabel('Real Signal');
a=sprintf('%s - Phase Sensitive Display',strrep(file, '_', '\_'));
title(a)

% plot imaginary part of signal
ax2=subplot(312);
plot(imag(fids(1:end)));
grid on;
xlabel('points');
ylabel('Imaginary Signal');

% plot magnitude of signal
ax3=subplot(313);
plot(abs(fids(1:end)));
grid on;
xlabel('points');
ylabel('Magnitude');

linkaxes([ax1,ax2,ax3],'x'); % link x axes