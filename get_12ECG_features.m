function features = get_12ECG_features(data, header_data)

       % addfunction path needed
        addpath(genpath('Tools/'))
        load('HRVparams_12ECG','HRVparams')
  p = gcp();
	% read number of leads, sample frequency and gain from the header.	

	[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

	HRVparams.Fs=Fs;
    HRVparams.PeakDetect.windows = floor(Total_time-1);
    HRVparams.windowlength = floor(Total_time);
    HRVparams.PVC.qrsth = 0.1;
	%parraler computing
  
    %prelocate
    ECG_ST_Elevations = NaN(1,12);
    ECG_Periods.PR = NaN;
    ECG_Periods.QS = NaN;
    ECG_Periods.QR = NaN;
    ECG_Periods.PT = NaN;
    ECG_Periods.TP = NaN;
    ECG_Periods.RAPR = NaN;
    PVC = -1;
    features_GEH = NaN(1,22);      
    AF_param = NaN(1,14);
    R_elevation = NaN(1,12);
    QRS_areas_vector = NaN(1,84);
    
    parfor i =1:num_leads
            Lead12wGain(i,:) = data(i,:)* gain(i);
    end
    
    %create high pass filter
    freq = 1;
    [b, a] = butter(6, freq / (Fs * 2), 'high');

    % median filter to remove bw
    parfor i=1:num_leads
        ECG_filt = filter(b,a,Lead12wGain(i,:)');
        ECG12filt(i,:) = medianfilter(ECG_filt, Fs);
    end
    
    try
        % convert 12Leads to XYZ leads using Kors transformation
        XYZLeads = Kors_git(ECG12filt);
        VecMag = vecnorm(XYZLeads');


        % Convert ECG waveform in rr intervals
        [~, rr, jqrs_ann, ~ , ~] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
        %sqi = [tSQI', SQIvalue'];

        % Find fiducial points using ECGKit
        ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);

    catch
        
    end

    f_pvc = parfeval(p,@get_PVC_feature,1,ECG12filt(7,:),HRVparams);
    
    
    if (exist('rr','var'))
        if (~isempty(rr))
            f_AF = parfeval(p,@get_AF_Parameters,1,rr,Fs);
        end
    end

    if (exist('Fid_pts','var'))
        f_GEH = parfeval(p,@get_GEH_feature,1,XYZLeads,Fid_pts,Fs);
        

        ECG_Periods = extract_ECG_Periods (Fid_pts,Fs);
        result = extract_ECG_ST_elevation_parameters (ECG12filt,Fid_pts);
        ECG_ST_Elevations = result.ST_elevation;
        R_elevation = result.R_elevation;
        QRS_areas_vector = compute_QRS_area (ECG12filt,Fid_pts,Fs);
    end

    try
        if (exist('f_pvc','var'))  
            PVC=fetchOutputs(f_pvc);
        end
        if (exist('f_AF','var')) 
            AF_param = fetchOutputs(f_AF);
        end
        if (exist('f_GEH','var'))
            features_GEH = fetchOutputs(f_GEH);
        end
        
    catch
    end
    ECG_Periods_vector = struct2vector (ECG_Periods);
    AF_important = [AF_param(2:7) AF_param(12) AF_param(14)];
    features = [age sex features_GEH AF_important PVC ECG_Periods_vector ECG_ST_Elevations QRS_areas_vector R_elevation];
            
end

function PVC = get_PVC_feature (ECG,HRVparams)
    try 
        PVC_flags = 0;
        QRS_numb = 0;
        
        [pvc_outputs,QRS] = PVC_detect(ECG,HRVparams);
        PVC_flags = PVC_flags + sum (pvc_outputs);
        QRS_numb = QRS_numb + length(QRS);
        
        PVC = PVC_flags / QRS_numb;
    catch
		PVC = -1;
    end 
end

function features = get_GEH_feature (XYZLeads,Fid_pts,Fs)
    features = NaN(1,22);
    try
        
        [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);
        GEH_features = GEH_analysis_git(XYZ_Median,Fid_pts_Median,Fs);

        features=GEH_features;

    catch

    end
end

function AF_param = get_AF_Parameters(rr,Fs)
    try    
        if (~isempty(rr)) 
            while (length (rr)<12)          
               rr = [rr rr]; 
            end
        if (length (rr)>59)
           rr = rr(1:59); 
        end

            AF_param = AF_features(round(rr*Fs),Fs);
        end
    catch
        AF_param = zeros(1,14);
    end
end

function ECG_Periods = extract_ECG_Periods (Fid_pts,Fs)
    Q = Fid_pts.Q;
    R = Fid_pts.R;
    S = Fid_pts.R;
    P = Fid_pts.P;
    T = Fid_pts.T;
    %If possible, compute periods (values in seconds).
    ECG_Periods.PR = nan (length (R),1);
    ECG_Periods.QS = nan (length (R),1);
    ECG_Periods.QR = nan (length (R),1);
    ECG_Periods.PT = nan (length (R),1);
    ECG_Periods.TP = nan (length (R),1);
    ECG_Periods.RAPR = nan (length (R)-1,1);
    
    for i=1:length (R)
        if ~(isnan(R(i))||isnan(P(i)))
            ECG_Periods.PR(i) = (R(i) - P(i))*(1/Fs);   
        end 
        
        if ~(isnan(S(i))||isnan(Q(i)))
            ECG_Periods.QS(i) = (S(i) - Q(i))*(1/Fs);   
        end
        
        if ~(isnan(R(i))||isnan(Q(i)))
            ECG_Periods.QR(i) = (R(i) - Q(i))*(1/Fs);   
        end
        
        if ~(isnan(P(i))||isnan(T(i)))
            ECG_Periods.PT(i) = (T(i) - P(i))*(1/Fs);   
        end
        
        if (i <length(R))%TP and RAPR use P_time from next cardiac cycle 
            if ~(isnan(P(i+1))||isnan(T(i)))
                ECG_Periods.TP(i) = (P(i+1) - T(i))*(1/Fs);   
            end
            if ~(isnan(R(i))||isnan(R(i+1))||isnan(ECG_Periods.PR(i)))
                ECG_Periods.RAPR (i)= ECG_Periods.PR(i) / ((R(i+1)-R(1))*(1/Fs));
            end
        end        
    end
end

function QRS_areas_vector = compute_QRS_area (ECG,Fid_pts,FS)
    n = 12;  
    retNR = 7;
    QRS_areas_vector = nan(1,n*retNR);

    QRSoff = Fid_pts.QRSoff;
    QRSon = Fid_pts.QRSon;
    for k=1:n
        QRS_areas_in_lead = nan (1,length(QRSoff));
        for i = 1:length(QRSoff)
            if (~(isnan(QRSon(i))||isnan(QRSoff(i))))
                QRS_areas_in_lead(i) = sum (ECG(k,QRSon(i):QRSoff(i))*(1/FS));
            end
        end
        p = (k-1)*retNR+1;
        QRS_areas_vector (p) =  min (QRS_areas_in_lead);
        QRS_areas_vector (p+1) =  max (QRS_areas_in_lead);
        QRS_areas_vector (p+2) =  median (rmmissing(QRS_areas_in_lead));
        QRS_areas_vector (p+3) =  std (rmmissing(QRS_areas_in_lead));
        
        [minDiff, maxDiff, meanDiff] = interbeatDifference (QRS_areas_in_lead);
        QRS_areas_vector (p+4) = minDiff;
        QRS_areas_vector (p+5) = maxDiff;
        QRS_areas_vector (p+6) = meanDiff;
    end
    
end

function result = extract_ECG_ST_elevation_parameters (ECG,Fid_pts)
    ECG_leads = 1:12;
    R_elevation = nan(1,length(ECG_leads));
    ST_elevation = zeros(1,length(ECG_leads));

    Toff = Fid_pts.Toff;    
    QRSoff = Fid_pts.QRSoff;
    QRSon = Fid_pts.QRSon;
    R = Fid_pts.R;
    
    for k=1:length(ECG_leads)
        nr = ECG_leads(k);
        
        d = (diff(ECG(nr,:)));
        QRSend =QRSoff;  
        QRSstart = QRSon;
        
        %fit QRS end to J point
        for a = 1:length(QRSend)
            if (~(isnan(QRSend(a))))
                if (ECG(nr,QRSend(a))< mean(ECG(nr,:))-(std(ECG(nr,:)*0.5)))
                    bp =0;
                    flag = true;
                    if (d(QRSend(a))>0)
                        znak = 1;
                    else
                        znak = 0;
                    end
                    while (flag)
                        QRSend(a) = QRSend(a) + 1;
                        bp = bp + 1;
                        if (bp >50)
                            break;
                        end
                        if (d(QRSend(a))<0 && znak == 1)
                            break;
                        end
                        if (d(QRSend(a))>0 && znak == 0)
                            break;
                        end
                    end
                end
            end
        end
        
        for a = 1:length(QRSstart)
            if (~(isnan(QRSstart(a))))
                if (ECG(nr,QRSstart(a))> mean(ECG(nr,:))+std(ECG(nr,:)))
                    bp =0;
                    flag = true;
                    if (d(QRSstart(a))>0)
                        znak = 1;
                    else
                        znak = 0;
                    end
                    while (flag)
                        QRSstart(a) = QRSstart(a) -1;
                        bp = bp + 1;
                        if (bp >50)
                            break;
                        end
                        if (d(QRSstart(a))<0 && znak == 1)
                            break;
                        end
                        if (d(QRSstart(a))>0 && znak == 0)
                            break;
                        end
                    end
                end
            end
        end
        J_Value = [];
        for lQ = 1: length (QRSend)
            if (isnan(QRSend(lQ)))
                J_Value(end+1) = nan;
            else
                J_Value(end+1) = (ECG(ECG_leads(k),(QRSend(lQ))));   
            end    
        end

        ECG_IzoLine = [nan];

        for j=2:length(Toff)
            if (~(isnan(Toff(j-1)) || isnan(QRSstart(j))))
                ECG_IzoLine (end+1) = median(ECG(ECG_leads(k),Toff(j-1):QRSstart(j)));
            else
                ECG_IzoLine (end+1) = nan;
            end
        end
        
        pom = J_Value-ECG_IzoLine;
        ST_elevation(k) = median(rmmissing(pom));
        if (~(isnan(R)))
            R_values = ECG(k,R);
            R_elevations = R_values-ECG_IzoLine;
            R_elevation(k) = median (rmmissing(R_elevations));
        end
    end
    result.ST_elevation = ST_elevation;
    result.R_elevation = R_elevation;
end

function [minDiff, maxDiff, meanDiff] = interbeatDifference (parameter)
    minDiff = nan;
    maxDiff = nan;
    meanDiff = nan; 
    if (length(rmmissing(parameter))>1)
        difference = diff(parameter);
        if (rmmissing(difference)>0)
            meanParameter = mean(rmmissing(parameter));
            interDiff = (difference ./ meanParameter)*100;
            minDiff = min (rmmissing(interDiff));
            maxDiff = max (rmmissing(interDiff));
            meanDiff = mean (rmmissing(interDiff));
        end
    end
end

function resultVector = struct2vector (S)
    fields = fieldnames(S);
    retNR = 7;
    resultVector =  nan(1,length(fields)*retNR);

    p = 1;
    for i=1:length(fields)
      value = S.(fields{i});
      if (~isempty(value))
          try
              resultVector (p) =  min (value);
              p = p + 1;
              resultVector (p) =  max (value);
              p = p + 1;
              resultVector (p) =  median (rmmissing(value));
              p = p + 1;
              resultVector (p) =  std (rmmissing(value));
              p = p + 1;
              [minDiff, maxDiff, meanDiff] = interbeatDifference (value);
              resultVector (p) =  minDiff;
              p = p + 1;
              resultVector (p) =  maxDiff;
              p = p + 1;
              resultVector (p) =  meanDiff;
              p = p + 1;
          catch
              p = (i+1)*retNR;
          end
      else
          p = p + retNR;
      end
    end
 
end
