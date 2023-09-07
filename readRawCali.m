function [cali_struct] = readRawCali(raw_path)
%Read in twix file or mrd file.
%   readRawCali(X) reads in the twix file or mrd file of the 129Xe
%   calibration scan and returns a structure containing all the variables
%   needed for processing.
% The data structure file contains:
%
% cali_struct.seq_name = sequence name (optional)
% cali_struct.weight = patient weight in pounds %
% cali_struct.te = TE in u-seconds;
% cali_struct.tr = TR in u-seconds;
% cali_struct.dwell_time = dwell time in nanoseconds;
% cali_struct.freq = gas excitation frequency in Hz;
% cali_struct.xeFreqMHz = gas excitation frequency in MHz;
% cali_struct.data = the FID data;
% cali_struct.scan_date = scan date in format (YYYY-MM-DD) optional;
% cali_struct.vref = reference voltage (V) optional;
% cali_struct.rf_excitation_ppm = rf excitation in ppm;

cali_struct = {};
[~, ~, file_extension] = fileparts(raw_path);

% Read in twix or P file and define associated variables
switch file_extension
    case '.dat'
        % Twix file from Siemens
        twix_obj = mapVBVD(raw_path);
        twix_obj.data = squeeze(double(twix_obj.image()));
        % If twix obj contains field sWipMemBlock, change to sWiPMemBlock (capital P)
        if isfield(twix_obj.hdr.MeasYaps, 'sWipMemBlock')
            temp = RenameField(twix_obj.hdr.MeasYaps, 'sWipMemBlock', 'sWiPMemBlock');
            twix_obj.hdr.MeasYaps = temp;
        end

        if isfield(twix_obj.hdr.Phoenix, 'sWipMemBlock')
            temp = RenameField(twix_obj.hdr.Phoenix, 'sWipMemBlock', 'sWiPMemBlock');
            twix_obj.hdr.Phoenix = temp;
        end
        % Extract variables and store in cali_struct
        npts = size(twix_obj.data, 1); % Number of samples per FID
        dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1, 1} * 10^-9; % Receiver bandwidth (kHz); works for both calibration types
        tr = twix_obj.hdr.Config.TR(1) * 1E-6; % Time between each sample
        fids = twix_obj.data;
        % Get the excitation frequency
        if isfield(twix_obj.hdr.Config, 'Frequency')
            % UVA Siemens File
            xeFreqMHz = twix_obj.hdr.Config.Frequency * 1e-6; %34.093484
        elseif isfield(twix_obj.hdr.Meas, 'lFrequency')
            % Duke Siemens File
            xeFreqMHz = twix_obj.hdr.Meas.lFrequency * 1e-6; %34.091516
        end
        % Get the scan date
        scanDate = twix_obj.hdr.Phoenix.tReferenceImage0;
        scanDate = strsplit(scanDate, '.');
        scanDate = scanDate{end};
        scanDateStr = [scanDate(1:4), '-', scanDate(5:6), '-', scanDate(7:8)];
        if isfield(twix_obj.hdr.Phoenix, 'sWiPMemBlock')
            % Duke twix file
            if isfield(twix_obj.hdr.Phoenix.sWiPMemBlock, 'adFree')
                VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{4};
                % seems to be in all sequences
                %rf_amp1 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;
                % in cali this is the 3rd flip angle (calibration, usually)
                %rf_amp2 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{3}.flAmplitude;
                %rf_amp3=1; %dummy override for when using hard pulse
            elseif isfield(twix_obj.hdr.Phoenix.sWiPMemBlock, 'alFree')
                % Duke twix file using UVA sequence
                VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1};
                % dissolved phase
                %rf_amp1 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{2}.flAmplitude;
                % calibration
                %rf_amp2 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{3}.flAmplitude;
            else
                disp('WARNING: twix file type not supported, cannot determine reference voltage')
            end
        else
            disp('WARNING: twix file type not supported, cannot determine reference voltage')
        end
        % Read RF excitation frequency
        mag_fstrength = twix_obj.hdr.Dicom.flMagneticFieldStrength; % Magnetic Field Strength
        excitation = twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1, 5}; % Read excitation from twix header,Cali version
        
        % RF excitation will be in ppm, likeley either 218ppm or 208 ppm at Duke
        gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
        rf_excitation_ppm = round(excitation/(gyro_ratio * mag_fstrength));

        % Save to output struct
        cali_struct.seq_name = twix_obj.hdr.Config.SequenceDescription;
        cali_struct.weight = twix_obj.hdr.Dicom.flUsedPatientWeight; %
        cali_struct.te = twix_obj.hdr.Phoenix.alTE{1};
        cali_struct.tr = twix_obj.hdr.Config.TR(1);
        cali_struct.dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1, 1}* 1E-9;
        cali_struct.freq = xeFreqMHz * 1e6;
        cali_struct.xeFreqMHz = xeFreqMHz;
        cali_struct.data = twix_obj.data;
        cali_struct.scan_date = scanDateStr;
        cali_struct.vref = VRef;
        cali_struct.rf_excitation_ppm = rf_excitation_ppm;

    case '.h5'
        % mrd file
        % Reading the k-space data
        dataset_data = h5read(raw_path, '/dataset/data');
        all_kspace_data = dataset_data.data;

        %Reshaping the k-space
        npts = size(all_kspace_data{1}, 1) / 2; % 512 for duke
        nfids = size(all_kspace_data, 1); % 520 for duke

        fids = [];
        for ii = 1:nfids
            fid = all_kspace_data{ii};
            i = 1;
            for j = 1:npts
                fids(j, ii) = double(complex(fid(i), fid(i+1)));
                i = i + 2;
            end
        end
        % Reading dwell_time from a k-space fid
        cali_data_head = dataset_data.head;
        dwell_time_all = cali_data_head.sample_time_us;
        dwell_time = double(dwell_time_all(1)) * 10^-6;

        % Dataset header variables are in the xml field - gives string
        dataset_header = h5read(raw_path, '/dataset/xml');

        xml_struct = xml2struct(dataset_header); % the xml2struct converts the string
        tr_mrd = str2double(xml_struct.ismrmrdHeader.sequenceParameters.TR.Text);

        % Threshold to not scale by 1E-3 some tr_mrd are in seconds
        if tr_mrd < 1E-1
            tr = tr_mrd; % this is 0.015 for duke
        else
            tr = tr_mrd * 1E-3;
        end

        vendor = xml_struct.ismrmrdHeader.acquisitionSystemInformation.systemVendor.Text;
        vendor = lower(vendor);
        switch vendor
            case 'ge'

            % conjugating the data solved the issues but we don't exactly know
            % the reason

            fids = conj(fids);
            % remove the first two points due to initial recovery time
            % and zero pad the end
            fids = fids(3:end, :);
            fids(size(fids, 1):size(fids, 1)+2, :) = 0;

            fids = fids / max(abs(fids(:)), [], 'all'); % max Normalizing
            xeFreqMHz_mrd = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 2}.value.Text);
            xeFreqMHz = xeFreqMHz_mrd * 1E-6; % 34.0923 MHz for duke
            % see if the excitation is in the mrd file
            try
                excitation = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
                % RF excitation will be in ppm, either 218ppm or 208 ppm
                rf_excitation = round(abs(excitation-xeFreqMHz*1E6)/(xeFreqMHz), 1);
            catch
                rf_excitation = 218;
            end
            % Philips - Cincinnati Children's Hospital Medical Center Data
            case 'philips'
            xeFreqMHz = str2double(xml_struct.ismrmrdHeader.userParameters.userParameterLong.value.Text) * 1E-6;
            rf_excitation = str2double(xml_struct.ismrmrdHeader.userParameters.userParameterDouble.value.Text);

            case 'siemens'
            xeFreqMHz_mrd = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 2}.value.Text);
            xeFreqMHz = xeFreqMHz_mrd * 1E-6;
            dwell_time = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1}.value.Text);
            % see if the excitation is in the mrd file
            try
                excitation = str2double(xml_struct.ismrmrdHeader.encoding.trajectoryDescription.userParameterDouble{1, 3}.value.Text);
                % RF excitation will be in ppm, either 218ppm or 208 ppm
                rf_excitation = round(abs(excitation-xeFreqMHz*1E6)/(xeFreqMHz), 1);
            catch
                rf_excitation = 218;
            end
            otherwise
            error('Unknown MR vendor')
        end
        
        
        if isfield(xml_struct.ismrmrdHeader, 'subjectInformation')
           if isfield(xml_struct.ismrmrdHeader.subjectInformation, 'patientWeightu_kg')
              cali_struct.weight = str2double(xml_struct.ismrmrdHeader.subjectInformation.patientWeightu_kg.Text);
           else
              cali_struct.weight = nan;
           end   
           if isfield(xml_struct.ismrmrdHeader.subjectInformation, 'patientID')
              cali_struct.seq_name = xml_struct.ismrmrdHeader.subjectInformation.patientID.Text;
           else
              cali_struct.seq_name = nan;
           end   
        else
           cali_struct.weight = nan;
           cali_struct.seq_name = nan;
        end   
        if strcmpi(vendor,'philips')
           cali_struct.te = str2double(xml_struct.ismrmrdHeader.sequenceParameters.TE.Text)*1e3;
        else
           cali_struct.te = str2double(xml_struct.ismrmrdHeader.sequenceParameters.TE.Text)*1e6;
        end
        cali_struct.tr = tr * 1e6;
        cali_struct.dwell_time = dwell_time;
        cali_struct.freq = xeFreqMHz * 1e6;
        cali_struct.xeFreqMHz = xeFreqMHz;
        cali_struct.data = fids;
        cali_struct.scan_date = nan;
        cali_struct.vref = nan;
        cali_struct.rf_excitation_ppm = rf_excitation;

    otherwise
        error('Unknown Raw File Type');
end
end
