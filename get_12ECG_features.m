function features = get_12ECG_features(data, header_data)

       % addfunction path needed
        addpath(genpath('Tools/'))
        load('HRVparams_12ECG','HRVparams')

	% read number of leads, sample frequency and gain from the header.	

	[recording,Total_time,num_leads,Fs,gain,age,sex]=extract_data_from_header(header_data);

	HRVparams.Fs=Fs;
        HRVparams.PeakDetect.windows = floor(Total_time-1);
        HRVparams.windowlength = floor(Total_time);
        HRVparams.PVC.qrsth = 0.1;
	%parraler computing
    p = gcp();

	ECG_PeriodsAndPeaks.ST_elevation = NaN;
    ECG_PeriodsAndPeaks.ECG_Periods = NaN;
    ECG_PeriodsAndPeaks.ST_data = NaN;
    features_GEH = NaN(1,24);
    features_GEH(1)=age;
    features_GEH(2)=sex;      
    AF_param = zeros(1,14);
        
        parfor i =1:num_leads
                Lead12wGain(i,:) = data(i,:)* gain(i);
        end


        % median filter to remove bw
        parfor i=1:num_leads
                ECG12filt(i,:) = medianfilter(Lead12wGain(i,:)', Fs);
        end
    try

        % convert 12Leads to XYZ leads using Kors transformation
        XYZLeads = Kors_git(ECG12filt);

        VecMag = vecnorm(XYZLeads');


        % Convert ECG waveform in rr intervals
        [t, rr, jqrs_ann, SQIvalue , tSQI] = ConvertRawDataToRRIntervals(VecMag, HRVparams, recording);
        sqi = [tSQI', SQIvalue'];

        % Find fiducial points using ECGKit
        ECG_header.nsig = 1; ECG_header.freq = Fs; ECG_header.nsamp = length(VecMag);
        wavedet_config.setup.wavedet.QRS_detection_only = 0;
        [Fid_pts,~,~] = wavedet_3D_ECGKit(VecMag', jqrs_ann', ECG_header, wavedet_config);

    catch
        ECG_PeriodsAndPeaks.ST_elevation = NaN;
        ECG_PeriodsAndPeaks.ECG_Periods = NaN;
        ECG_PeriodsAndPeaks.ST_data = NaN;
		features_GEH = NaN(1,24);
        features_GEH(1)=age;
        features_GEH(2)=sex;      
        AF_param = zeros(1,14);
	end

    
    if (exist('rr','var'))
        if (~isempty(rr))
            f_AF = parfeval(p,@get_AF_Parameters,1,rr,Fs);
        end
    end

    if (exist('Fid_pts','var'))
        f_st = parfeval(p,@extract_ECG_Periods,1,ECG12filt,Fid_pts,Fs);
        f_GEH = parfeval(p,@get_GEH_feature,1,XYZLeads,Fid_pts,Fs, age, sex);
    end
    try 
        PVC = get_PVC_feature (ECG12filt(7,:),HRVparams);

        if (exist('f_st','var'))   
           ECG_PeriodsAndPeaks = fetchOutputs(f_st);
        end
        if (exist('f_GEH','var')) 
            features_GEH = fetchOutputs(f_GEH);
        end
        if (exist('f_AF','var')) 
            AF_param = fetchOutputs(f_AF);
        end
    catch
       PVC = -1; 
    end
    
    result.GEH = features_GEH;
    result.AF_param = AF_param;
    result.PVC = PVC;
    result.resultVector = vectorization (ECG_PeriodsAndPeaks.ST_data,ECG_PeriodsAndPeaks.ECG_Periods);

    features = [features_GEH AF_param PVC  result.resultVector];
            
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

function features = get_GEH_feature (XYZLeads,Fid_pts,Fs, age, sex)
    try
        [XYZ_Median,Fid_pts_Median] = Time_coherent_code_github(XYZLeads,Fid_pts,Fs);
        GEH_features = GEH_analysis_git(XYZ_Median,Fid_pts_Median,Fs);
        features(1)=age;
        features(2)=sex;
        features(3:24)=GEH_features;
    catch
        features = NaN(1,24);
        features(1)=age;
        features(2)=sex;
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

function result = extract_ECG_Periods (ECG,Fid_pts,Fs)
    ECG_leads = 1:12;
    ST_elevation = zeros(length(ECG_leads),1);
    R_sign = zeros(length(ECG_leads),1);
    ST_data={};
    
    
    QRSoff = Fid_pts.QRSoff;
    QRSon = Fid_pts.QRSon;
    Toff = Fid_pts.Toff;
    R = Fid_pts.R;
    P = Fid_pts.P;
    
    %If possible, compute period (values are in seconds).
    ECG_Periods.PR = zeros (length (R),1);
    ECG_Periods.QS = zeros (length (R),1);
    ECG_Periods.QR = zeros (length (R),1);
    for i=1:length (R)
        if (isnan(R(i))||isnan(P(i)))
            ECG_Periods.PR(i) = nan;
        else
            ECG_Periods.PR(i) = (R(i) - P(i))*(1/Fs);   
        end 
        
        if (isnan(QRSoff(i))||isnan(QRSon(i)))
            ECG_Periods.QS(i) = nan;
        else
            ECG_Periods.QS(i) = (QRSoff(i) - QRSon(i))*(1/Fs);   
        end
        
        if (isnan(R(i))||isnan(QRSon(i)))
            ECG_Periods.QR(i) = nan;
        else
            ECG_Periods.QR(i) = (R(i) - QRSon(i))*(1/Fs);   
        end
    end
    
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

        R_peak_value = mean(ECG(nr,rmmissing(R)))-mean(rmmissing(ECG_IzoLine));
        if (R_peak_value>0)
            R_sign (k) = 1;
        else
            R_sign (k) = -1;
        end

        pom = J_Value-ECG_IzoLine;
        ST_elevation(k) = median(rmmissing(pom));
        
        R_peak = zeros (1,length(R));
        for b=1:length(R)
            if isnan(R(b))
                R_peak(b) = nan;
            else
                R_peak(b) = ECG(nr,R(b));
            end
        end
        element = struct('ST_elevations',pom, 'J_Value',J_Value, 'ECG_IzoLine',ECG_IzoLine, 'R_peak',R_peak);
 
        
        ST_data{k} = element;
    end
       
    result.ST_elevation = ST_elevation.*R_sign;
    result.ST_data = ST_data;
    result.ECG_Periods = ECG_Periods;
end

function resultVector = vectorization (ECG_Values,ECG_Periods)
    resultVector = zeros(1,108);
    p = 1;
    if (length(ECG_Values)==12)
        for n=1:12
              ST_elevations = rmmissing(ECG_Values{n}.ST_elevations); 
              ECG_IzoLine = rmmissing(ECG_Values{n}.ECG_IzoLine); 
              if (~isempty(ST_elevations))
                  resultVector (p) =  min (ST_elevations);
                  p = p + 1;
                  resultVector (p) =  max (ST_elevations);
                  p = p + 1;
                  resultVector (p) =  median (ST_elevations);
                  p = p + 1;
                  resultVector (p) =  std (ST_elevations);
                   p = p + 1;
              else
                   resultVector (p) = nan;
                   p = p + 1;
                   resultVector (p) = nan;
                   p = p + 1;
                   resultVector (p) = nan;
                   p = p + 1;
                   resultVector (p) = nan;
                   p = p + 1;
              end

              if (~isempty(ECG_IzoLine))
                  resultVector (p) =  min (ECG_IzoLine);
                  p = p + 1;
                  resultVector (p) =  max (ECG_IzoLine);
                  p = p + 1;
                  resultVector (p) =  median (ECG_IzoLine);
                  p = p + 1;
                  resultVector (p) =  std (ECG_IzoLine);
                  p = p + 1;
              else
                   resultVector (p) = nan;
                   p = p + 1;
                   resultVector (p) = nan;
                   p = p + 1;
                   resultVector (p) = nan;
                   p = p + 1;
                   resultVector (p) = nan;
                   p = p + 1;
              end
        end
    else
       for i=1:96
           resultVector (p) = nan;
           p = p + 1;
       end
    end
      if (isstruct(ECG_Periods))
          PR = rmmissing(ECG_Periods.PR);
          if (~isempty(PR))
              resultVector (p) =  min (PR);
              p = p + 1;
              resultVector (p) =  max (PR);
              p = p + 1;
              resultVector (p) =  median (PR);
              p = p + 1;
              resultVector (p) =  std (PR);
              p = p + 1;
          else
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
          end

          QS = rmmissing(ECG_Periods.QS);
          if (~isempty(QS))
              resultVector (p) =  min (QS);
              p = p + 1;
              resultVector (p) =  max (QS);
              p = p + 1;
              resultVector (p) =  median (QS);
              p = p + 1;
              resultVector (p) =  std (QS);
              p = p + 1;
          else
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
          end
          QR = rmmissing(ECG_Periods.QR);
          if (~isempty(QR))
              resultVector (p) =  min (QR);
              p = p + 1;
              resultVector (p) =  max (QR);
              p = p + 1;
              resultVector (p) =  median (QR);
              p = p + 1;
              resultVector (p) =  std (QR);
              p = p + 1;
          else
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
               resultVector (p) = nan;
               p = p + 1;
          end
      else
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;
           p = p + 1;
           resultVector (p) = nan;            
      end
end

