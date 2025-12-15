function z = zHeartbeat(o,filename,samplingRate)
    % Read heartbeat from a recording in the MIT-BIH dataset 
    % https://physionet.org/content/mitdb/1.0.0/
    % and return the impedance to use at eacht timepoint of an ArtSim
    % simulation.
    %
    % The example data set (MIT_BIH_100.csv) has two leads,
    % was recorded at 360 samples/s from a person with a heartbeat of ~ 75 bpm.
    % This function uses the MLII lead.
    %
    % This function models the balistocardiac effect as a delayed, damped 
    % sinusoid induced by the ventricular contraction (R-wave in the ECG). 
    % 
    % The output (a factor to scale the impedance) has a mean of 1 and a
    % maximum of o.zAmplitude
    %
    % By specifying zFun = @zHeartbeat in an ArtSim object, this function 
    % is used to simulate impedance artifacts.
    % 
    arguments
        o (1,1) artSim
        filename (1,1) string = "MIT_BIH_100.csv" 
        samplingRate =360;
    end
    
    assert(exist(filename,"file"),"MIT BIH file not found: %s",filename);
    warning('off','MATLAB:table:ModifiedAndSavedVarnames')
    ECG = readtable(filename);
    ecgTime  =seconds(ECG.x_sample__/samplingRate);    
    warning('on','MATLAB:table:ModifiedAndSavedVarnames');   
    assert(seconds(o.tSimulate(end))<ecgTime(end),"Simulation is too long for this heartbeat recording")
    % Convert samples to seconds
    ecg = ECG.x_MLII_ -median(ECG.x_MLII_); 
    ecg(ecg<0) = 0; % Removes smaller (P,T) waves and keeps the ventricular contraction (R-peak).

    % Convolve the R-peaks with a damped sinusoid to model the ballistic
    % artifact
    bodyResonanceFrequency = 5;      % Wobble frequency [Hz]
    tau = 0.1;                       % Time constant of the dampening [s]
    delay = 0.2;                    % Wobble follows the ECG peak by this  [s]        
    t= (0:1/samplingRate:1);
    nrDelaySamples = round(delay*samplingRate);
    ballisticModel  = [zeros(1,nrDelaySamples) exp(-t/tau).*sin(2*pi*bodyResonanceFrequency*t)];    
    b =  conv(ecg,ballisticModel,"full"); % Convolution model

    B = timetable(ecgTime,b(1:numel(ecg)),'VariableNames',"ballistic");
    % Resample to match the simulation 
    B = retime(B,seconds(o.tSimulate),"spline");
    z = B.ballistic- mean(B.ballistic);
    z = z./max(z);% Scale <=1
    z= 1+o.zAmplitude*z;  % Scale to zAmplitude
end