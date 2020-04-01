function [features, AF_param, PVC] = get_12ECG_features(data, header_data)

       % addfunction path needed
        addpath(genpath('Tools/'))
        load('HRVparams_12ECG','HRVparams')

	% read number of leads, sample frequency and gain from the header.	

	[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

	HRVparams.Fs=Fs;
        HRVparams.PeakDetect.windows = floor(Total_time-1);
        HRVparams.windowlength = floor(Total_time);
    HRVparams.PVC.qrsth = 0.1;
	

    parfor i =1:num_leads
            Lead12wGain(i,:) = data(i,:)* gain(i);
    end

    
    % median filter to remove bw
    parfor i=1:num_leads
            ECG12filt(i,:) = medianfilter(Lead12wGain(i,:)', Fs);
    end

    try 
        PVC_flags = 0;
        QRS_numb = 0;
        %select = [7, 8, 9, 10];
        %for i=1:length(select)
            [pvc_outputs,QRS] = PVC_detect(ECG12filt(7,:),HRVparams);
            PVC_flags = PVC_flags + sum (pvc_outputs);
            QRS_numb = QRS_numb + length(QRS);
        %end
        PVC = PVC_flags / QRS_numb;
    catch
		PVC = 0;
    end    
    try          
        % convert 12Leads to XYZ leads using Kors transformation
        XYZLeads = Kors_git(ECG12filt);

        VecMag = vecnorm(XYZLeads');
        
        % Convert ECG waveform in rr intervals
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
        sqi = [tSQI', SQIvalue'];
        
        if (~isempty(rr)) 
            while (length (rr)<12)          
               rr = [rr rr]; 
            end
            if (length (rr)>59)
               rr = rr(1:59); 
            end

            AF_param = AF_features(round(rr*Fs),Fs);
        end


        %
        % Find fiducial points using ECGKit
        ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);

        [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);

        GEH_features = GEH_analysis_git(XYZ_Median,Fid_pts_Median,Fs);

        features(1)=age;
        features(2)=sex;
        features(3:24)=GEH_features;


    catch
		features = NaN(1,24);
        AF_param = zeros(1,14);
    end

end

