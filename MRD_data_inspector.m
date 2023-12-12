% Extract key header variables and plot raw data from MRD files.

%% Set flag arguments
clear; close all;              % clear current variables and close figures
first_frames=5;                % number of initial fids to plot typically 5-20

%% Select MRD file
[file, path] = uigetfile('*.*', 'Select file'); % open UI for selecting file, starts in current dir
file_with_path = strcat(path, file); % join path and file to open
dataset = ismrmrd.Dataset(file_with_path, 'dataset');
ismrmrd_header = ismrmrd.xml.deserialize(dataset.readxml);
%% Extract key header variables

% convert user parameter fields to maps for easy query

if contains(file, "dixon") || contains(file, "calibration")
    general_user_params_long = containers.Map();
    data_struct = ismrmrd_header.userParameters.userParameterLong;
    for i = 1:numel(data_struct)
        general_user_params_long(data_struct(i).name) = data_struct(i).value;
    end
end

if contains(file, "proton") || contains(file, "dixon")
    general_user_params_string = containers.Map();
    data_struct = ismrmrd_header.userParameters.userParameterString;
    for i = 1:numel(data_struct)
        general_user_params_string(data_struct(i).name) = data_struct(i).value;
    end

    traj_description_user_params = containers.Map();
    data_struct = ismrmrd_header.encoding.trajectoryDescription.userParameterLong;
    for i = 1:numel(data_struct)
        traj_description_user_params(data_struct(i).name) = data_struct(i).value;
    end
end

% general variables
study_date = ismrmrd_header.studyInformation.studyDate; % in YYYY-MM-DD
try
    patient_id = ismrmrd_header.subjectInformation.patientID;
catch
    patient_id = nan;
    disp("No patient ID found.")
end
vendor = ismrmrd_header.acquisitionSystemInformation.systemVendor;
institution = ismrmrd_header.acquisitionSystemInformation.institutionName;
field_strength = ismrmrd_header.acquisitionSystemInformation.systemFieldStrength_T; % in T
te90 = ismrmrd_header.sequenceParameters.TE; % in us
sample_time = dataset.readAcquisition().head.sample_time_us(1); % in us
fov = [ismrmrd_header.encoding.reconSpace.fieldOfView_mm.x, ...
    ismrmrd_header.encoding.reconSpace.fieldOfView_mm.y, ...
    ismrmrd_header.encoding.reconSpace.fieldOfView_mm.z]; % in mm

% scan-specific variables
if contains(file, 'proton')
    tr_proton = ismrmrd_header.sequenceParameters.TR(1); % in ms
    flip_angle_proton = ismrmrd_header.sequenceParameters.flipAngle_deg(1); % in degrees
    matrix_size_z = ismrmrd_header.encoding.reconSpace.matrixSize.z;
    orientation = general_user_params_string("orientation");
    ramp_time = traj_description_user_params("ramp_time"); % in us
elseif contains(file, 'dixon')
    freq_center = general_user_params_long("xe_center_frequency"); % in Hz
    freq_dis_excitation_hz = general_user_params_long("xe_dissolved_offset_frequency"); % in Hz
    tr_gas = ismrmrd_header.sequenceParameters.TR(1);
    tr_dissolved = ismrmrd_header.sequenceParameters.TR(2);
    flip_angle_gas = ismrmrd_header.sequenceParameters.flipAngle_deg(1);
    flip_angle_dissolved = ismrmrd_header.sequenceParameters.flipAngle_deg(2);
    matrix_size_z = ismrmrd_header.encoding.reconSpace.matrixSize.z;
    orientation = general_user_params_string("orientation");
    ramp_time = traj_description_user_params("ramp_time"); % in us
elseif contains(file, "calibration")
    freq_center = general_user_params_long("xe_center_frequency"); % in Hz
    freq_dis_excitation_hz = general_user_params_long("xe_dissolved_offset_frequency"); % in Hz
    tr_gas = ismrmrd_header.sequenceParameters.TR(1);
    tr_dissolved = ismrmrd_header.sequenceParameters.TR(2);
    flip_angle_gas = ismrmrd_header.sequenceParameters.flipAngle_deg(1);
    flip_angle_dissolved = ismrmrd_header.sequenceParameters.flipAngle_deg(2);
end

% calculate RF excitation in ppm
if contains(file, "dixon") || contains(file, "calibration")
    gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
    rf_excitation = freq_dis_excitation_hz/(gyro_ratio * field_strength);
end

%% Extract FID acquisition data and trajectories

% read k-space data
% fids = cell2mat(dataset.readAcquisition().data);
npts = size(dataset.readAcquisition(1).data{1},1);
nfids = dataset.getNumberOfAcquisitions;
fids_cell = dataset.readAcquisition().data;
fids = zeros(nfids, npts);
for i=1:nfids
    fids(i,:) = transpose(double(fids_cell{i}));
end

% if data from GE scanner, take complex conjugate
if strcmpi(vendor,'ge')
    fids = conj(fids);
end

% read trajectory data
if contains(file, "proton") || contains(file, "dixon")
    traj_cell = dataset.readAcquisition().traj;
    traj = zeros(nfids, npts, 3);
    for i=1:nfids
        traj(i,:,:) = transpose(double(traj_cell{i}));
    end
end

% get contrast and bonus spectra labels
contrast_labels = dataset.readAcquisition().head.idx.contrast;
bonus_spectra_labels = dataset.readAcquisition().head.measurement_uid;

%% print out key variables

% general variables
fprintf('\tFile = %s\n',file);
fprintf('\tPatient ID = %s\n',patient_id);
fprintf('\tStudy date = %s\n',study_date);
fprintf('\tSystem vendor = %s\n',vendor);
fprintf('\tInstitution = %s\n',institution);
fprintf('\tField strength = %0.2f T\n',field_strength);
fprintf('\tTE90 = %0.3f ms\n',te90);
fprintf('\tSample time = %0.2f us\n',sample_time);
fprintf('\tFOV (x) = %0.0f mm\n', fov(1));
fprintf('\tFOV (y) = %0.0f mm\n', fov(2));
fprintf('\tFOV (z) = %0.0f mm\n', fov(3));
fprintf('\tNum FIDs = %0.0f\n',nfids);
fprintf('\tNum pts in FID (npts) = %0.0f\n',npts);
fprintf('\tNum bonus spectra = %0.0f\n',sum(bonus_spectra_labels==1));

% scan-specific variables
if contains(file, 'proton')
    fprintf('\tTR (proton) = %0.2f ms\n',tr_proton);
    fprintf('\tFlip angle (proton) = %0.0f deg\n',flip_angle_proton);
    fprintf('\tMatrix size (z) = %0.0f\n', matrix_size_z);
    fprintf('\tOrientation = %s\n',orientation);
    fprintf('\tRamp time = %0.0f us\n',ramp_time);
elseif contains(file, 'dixon')
    fprintf('\tCenter frequency = %0.0f Hz\n',freq_center);
    fprintf('\tOffset frequency = %0.0f Hz\n',freq_dis_excitation_hz);
    fprintf('\tRF excitation = %0.1f ppm\n',rf_excitation);
    fprintf('\tTR (gas) = %0.1f ms\n',tr_gas);
    fprintf('\tTR (dissolved) = %0.1f ms\n',tr_dissolved);
    fprintf('\tFlip angle (gas) = %0.2f deg\n',flip_angle_gas);
    fprintf('\tFlip angle (dissolved) = %0.0f deg\n',flip_angle_dissolved);
    fprintf('\tMatrix size (z) = %0.0f\n', matrix_size_z);
    fprintf('\tOrientation = %s\n',orientation);
    fprintf('\tRamp time = %0.0f us\n',ramp_time);
elseif contains(file, "calibration")
    fprintf('\tCenter frequency = %0.0f Hz\n',freq_center);
    fprintf('\tOffset frequency = %0.0f Hz\n',freq_dis_excitation_hz);
    fprintf('\tRF excitation = %0.1f ppm\n',rf_excitation);
    fprintf('\tTR (gas) = %0.1f ms\n',tr_gas);
    fprintf('\tTR (dissolved) = %0.1f ms\n',tr_dissolved);
    fprintf('\tFlip angle (gas) = %0.2f deg\n',flip_angle_gas);
    fprintf('\tFlip angle (dissolved) = %0.0f deg\n',flip_angle_dissolved);
end

%% Plot a select number of first FIDs sequentially

% check if number of requested fids exceeds number of total fids
if first_frames>nfids
    first_frames=nfids;
    sprintf('\nNumber of plot frames exceeded. Setting to %g\n',first_frames);
end 

% plot first FIDs
figure();
plot(abs(fids(1:first_frames*npts)));

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

% plot all FIDs sequentially
figure();
plot(abs(fids(1:end))); 

% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All FIDs Sequential', strrep(file, '_', '\_'));
title(a)
box on

%% Plot all FIDs overlapping (to look for noise)

% plot all FIDs overlapping
figure();
plot(abs(fids)); 

% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All Dixon FIDS Overlapping', strrep(file, '_', '\_'));
title(a)
box on

%% Plot all gas FIDs overlapping (minus bonus spectra)

% plot gas FIDs overlapping without bonus spectra
figure();
plot(abs(fids(:,contrast_labels==1 & bonus_spectra_labels==0)))

% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All Gas FIDS Overlapping (minus bonus spectra)', strrep(file, '_', '\_'));
title(a)
box on

%% Plot all dissolved FIDs overlapping (minus bonus spectra)

% plot gas FIDs overlapping without bonus spectra
figure();
plot(abs(fids(:,contrast_labels==2 & bonus_spectra_labels==0)))

% set title and axis labels
xlabel('points')
ylabel('magnitude')
a=sprintf('%s All Dissolved FIDS Overlapping (minus bonus spectra)', strrep(file, '_', '\_'));
title(a)
box on

%% Plot real and imaginary parts of signal

% initialize figure
figure();
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