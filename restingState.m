%[text] %[text:anchor:T_360F0929] # Artifacts in resting state EEG data
%[text] This script uses ArtSim's option to load data (from a Fieldtrip or EEGLAB data structure) and investigates whether (simulated) tACS can be removed and the true underlying signals recovered. 
%[text] Specifically, I used 4 minutes of (eyes open) resting state activity recorded at electrodes near Oz.  In these data the 1/f decay rate as well as the strength of the alpha rhythm varies considerably across participants. 
%[text] The simulated tACS signal applies a 1 mA current at 10 Hz  that results in a 15 mV artifact in the EEG signal.
%[text] I then apply the three artifiact removal algorithms and compute the residual artifact (expressed as the percent error in the estimated power) near 10 Hz (directly affected by tACS) and near 5 Hz (control frequency).
%[text] %[text:anchor:H_D32E8D0B] ## Setup Simulation
%[text] First define the stimulation and recording parameters.
startup; % Make sure all toolboxes are on the path

defaults.reproducible              = true;    % Fix the random number generator to get reproducible results
%[text] Stimulation and Recording Equipment:
defaults.stimulatorSamplingRate    = 1e3;     % 1 kHz 
defaults.currentResolution         = 1e-6;    % 1 muA
defaults.highPass                  = 0.01;    % Analog filter of the recording amplifier
defaults.lowPass                   = 100;
defaults.recordingSamplingRate     = 250;
%[text] Experiment:
defaults.duration                  = 4*60;   % 4 min of data
defaults.tacsAmplitude              = 1e-3;    % 1 mA 
defaults.resistance                 = 15;     % 15 mV artifact per mA current
defaults.tacsFrequency              = 10;      % 5 Hz tACS
%[text] And setup the base parameters for artifiact removal
arDefaults.slide                     = false;   % No sliding segment
arDefaults.slack                     = 5e-3;    % 5 ms to realign segments
arDefaults.nrSegsPerWindow           = 10;      % 10 segments to average artifact
arDefaults.pcaNr                     = 2;       % First two PC are removed
arDefaults.pcaVarExplained           = Inf;     % Ignore explained variance, always remove 2.
arDefaults.pcaNrSegsPerWindow        = Inf;     % Use all segments (plus the tCS signal; below) to determine basis set.
arDefaults.referenceSegment          = 1;       % Align to first segment
arDefaults.recordingSamplingRate     = defaults.recordingSamplingRate; % Match to the recording
arDefaults.segmentDuration           = 3;

%%
%[text] %[text:anchor:H_DDD99D7B] ## Estimate Alpha Power
%[text] Create the artSim object using the default parameters defined above.
A                   = artSim(defaults);
ar                          = arDefaults;
ar.tacsFrequency            = A.tacsFrequency;
arModes = ["RBAR" "FBAR" "FASTR"];


pctError = table();
if ~isempty(getenv("DJ_HOST")) && false %[output:group:6f8eda59]
    %% kLab only: retrieve the data from the datajoint server,
    allRestData = ns.Epoch & 'dimension="rest"' & (ns.Experiment & 'paradigm ="rvsRestingState"');    
    sCntr =0;
    for thisRestData = fetch(allRestData)'
        sCntr = sCntr+1;
        % Extract a Fieldtrip dataset
        ft = ephys.fieldtrip.dataset(thisRestData,dimension="rest");
        % Add the data to the ArtSim object (focus on electrodes near Oz)
        eeg = readData(A,ft,format="FIELDTRIP",channel= [125 126 138 139 117]);         
        % Save the data for future use without database access.
        save("../data/" + thisRestData.subject  + ".rest.mat","eeg");
        for arMode = arModes               
            thisPctError = cleanAndEstimatePower(A,ar,arMode,thisRestData.subject);
            pctError = [pctError;thisPctError]; %#ok<AGROW>
        end
    end
else
    % Read subject.rest.mat data from ../data created by the code abobe
    subjects = dir('../data/*.rest.mat');
    for sCntr = 1:numel(subjects)
        load(fullfile(subjects(sCntr).folder,subjects(sCntr).name),"eeg");
        % Add the data to the ArtSim object (focus on electrodes near Oz)
        readData(A,eeg,format="RAW",channel= [125 126 138 139 117]);  %[output:6e76908b] %[output:18e88981] %[output:9afbfc7d] %[output:2733656e] %[output:7c45ded7] %[output:4ecc6ffb] %[output:2843b0ba] %[output:544bbf6c] %[output:8f581fb7] %[output:6aae124b] %[output:4336304a] %[output:23985348] %[output:956c781d] %[output:0a3e506b] %[output:4efbe68f] %[output:734af950] %[output:0f211e6e] %[output:93482ee6] %[output:338c7722] %[output:565458d5] %[output:01fc80a2] %[output:5a393e9a] %[output:72713035] %[output:8927da9d] %[output:6ead7dc4] %[output:396b2734] %[output:46087d98]
        for arMode = arModes               
            thisPctError = cleanAndEstimatePower(A,ar,arMode,extractBefore(subjects(sCntr).name,".rest.mat")); %[output:1b4cbd2f] %[output:1b89b174] %[output:059aba78] %[output:8272f826] %[output:4da93917] %[output:90086f90] %[output:58abf7e3] %[output:2969cf7a] %[output:62b6e94d] %[output:4020842c] %[output:4d0341e7] %[output:3f99c3a5] %[output:8e19b20b] %[output:4fbbc53e] %[output:9d45243f] %[output:503287e7] %[output:8d9a8b7a] %[output:90f6b5f8] %[output:682e5da1] %[output:0f5f2451] %[output:69cb7e61] %[output:6ea0daa2] %[output:387be710] %[output:4efaf9c7] %[output:44da7429] %[output:9534d4b4] %[output:136ec0c3] %[output:63acea6d] %[output:35f204e6] %[output:27020378] %[output:3f45bc1b] %[output:128be23c] %[output:917a59cd] %[output:5b0aaa32] %[output:14085279] %[output:469d40ea] %[output:28892c8b] %[output:10c20afb] %[output:43df5a9d] %[output:4c23c1a8] %[output:611a36b2] %[output:4032f9d3] %[output:8fcde6f3] %[output:480a9f12] %[output:9cc25b1c] %[output:00562850] %[output:128cbb43] %[output:8736beca] %[output:77c78d20] %[output:99b43421] %[output:19e4bff7] %[output:973ec693] %[output:5923081a] %[output:44dc52b1] %[output:6afd1434] %[output:2e1b636a] %[output:48a354c4] %[output:5d8e6577] %[output:572d0ab0] %[output:8d42e345] %[output:4869747d] %[output:08b15ee4] %[output:237161ab] %[output:918abdb6] %[output:90728558] %[output:5125bbce] %[output:2fe8c8a1] %[output:5d2383e2] %[output:1c7d796a] %[output:4d18e78a] %[output:76589a7a] %[output:5180fd7c] %[output:2e974166] %[output:468d7a24] %[output:9da770da] %[output:83190844] %[output:75e301c5] %[output:4af4b559] %[output:3105f2b4] %[output:23f9252a] %[output:72f65666]
            pctError = [pctError;thisPctError]; %#ok<AGROW>
        end
    end
end %[output:group:6f8eda59]
%%
%[text] ## Summarize the results
figure %[output:70ab0732]
tiledlayout(3,2); %[output:70ab0732]
for arMode =arModes %[output:group:5f925c9f]
    for e = ["e5" "e10"]
        nexttile %[output:70ab0732]
            thisE = pctError{pctError.arMode==arMode,e};
            histogram(thisE,BinMethod="fd");
            title (arMode +"/" + e + ": " + num2str(median(thisE),2) + "% SD: " + num2str(std(thisE),2) + "%");
    end
end %[output:group:5f925c9f]
%%
%[text] ## Helper Functions
%[text] Simple example analysis; remove artifacts, then estimate the power at 10 Hz and 5 Hz, and return the errors as a percentage.
function pctError = cleanAndEstimatePower(o,ar,arMode,subject)

r  = simPlotEEG(o,ar,arMode,freqRes =1,axs=[]);
s10 = r.amplitude(r.frequency>9 & r.frequency <11,1:2);
s5 = r.amplitude(r.frequency>4 & r.frequency <6,1:2);            
thisPct5 = mean(100*abs((s5(:,2)-s5(:,1))./s5(:,1)),"omitmissing");
thisPct10 = mean(100*abs((s10(:,2)-s10(:,1))./s10(:,1)),"omitmissing");        
pctError = table(string(subject),arMode,thisPct5,thisPct10,'VariableNames',{'subject','arMode','e5','e10'});
end


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:6e76908b]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:1b4cbd2f]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.02 rmse = 66\/0 muV ae = 23.2%)\nAllocating memory...\n","truncated":false}}
%---
%[output:1b89b174]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:059aba78]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.02 rmse = 47\/0 muV ae = 123.3%)\n","truncated":false}}
%---
%[output:18e88981]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:8272f826]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=0.99 rmse = 0\/0 muV ae = 0.4%)\nFBAR (r=0.05 rmse = 66\/0 muV ae = 22.0%)\nAllocating memory...\n","truncated":false}}
%---
%[output:4da93917]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:90086f90]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.04 rmse = 47\/0 muV ae = 69.4%)\n","truncated":false}}
%---
%[output:9afbfc7d]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:58abf7e3]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=-0.01 rmse = 66\/0 muV ae = 31.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:2969cf7a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:62b6e94d]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.00 rmse = 47\/0 muV ae = 97.1%)\n","truncated":false}}
%---
%[output:2733656e]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:4020842c]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=-0.01 rmse = 66\/0 muV ae = 38.6%)\nAllocating memory...\n","truncated":false}}
%---
%[output:4d0341e7]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:3f99c3a5]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=-0.00 rmse = 47\/0 muV ae = 135.2%)\n","truncated":false}}
%---
%[output:7c45ded7]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:8e19b20b]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 1.1%)\nFBAR (r=0.83 rmse = 66\/0 muV ae = 9.1%)\nAllocating memory...\n","truncated":false}}
%---
%[output:4fbbc53e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:9d45243f]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.78 rmse = 84\/0 muV ae = 97.1%)\n","truncated":false}}
%---
%[output:4ecc6ffb]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:503287e7]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.4%)\nFBAR (r=0.08 rmse = 66\/0 muV ae = 22.1%)\nAllocating memory...\n","truncated":false}}
%---
%[output:8d9a8b7a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:90f6b5f8]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.06 rmse = 47\/0 muV ae = 125.1%)\n","truncated":false}}
%---
%[output:2843b0ba]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:682e5da1]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.76 rmse = 66\/0 muV ae = 15.6%)\nAllocating memory...\n","truncated":false}}
%---
%[output:0f5f2451]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:69cb7e61]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.77 rmse = 43\/0 muV ae = 34.7%)\n","truncated":false}}
%---
%[output:544bbf6c]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:6ea0daa2]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=0.99 rmse = 0\/0 muV ae = 0.6%)\nFBAR (r=0.21 rmse = 66\/0 muV ae = 18.0%)\nAllocating memory...\n","truncated":false}}
%---
%[output:387be710]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:4efaf9c7]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.22 rmse = 47\/0 muV ae = 54.7%)\n","truncated":false}}
%---
%[output:8f581fb7]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:44da7429]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.1%)\nFBAR (r=0.59 rmse = 66\/0 muV ae = 10.9%)\nAllocating memory...\n","truncated":false}}
%---
%[output:9534d4b4]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:136ec0c3]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.52 rmse = 39\/0 muV ae = 118.1%)\n","truncated":false}}
%---
%[output:6aae124b]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:63acea6d]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.16 rmse = 66\/0 muV ae = 19.2%)\nAllocating memory...\n","truncated":false}}
%---
%[output:35f204e6]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:27020378]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.15 rmse = 47\/0 muV ae = 82.4%)\n","truncated":false}}
%---
%[output:4336304a]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:3f45bc1b]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 8.1%)\nFBAR (r=0.96 rmse = 66\/0 muV ae = 31.2%)\nAllocating memory...\n","truncated":false}}
%---
%[output:128be23c]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:917a59cd]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.94 rmse = 50\/0 muV ae = 360.1%)\n","truncated":false}}
%---
%[output:23985348]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:5b0aaa32]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.06 rmse = 66\/0 muV ae = 13.1%)\nAllocating memory...\n","truncated":false}}
%---
%[output:14085279]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:469d40ea]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.06 rmse = 47\/0 muV ae = 82.6%)\n","truncated":false}}
%---
%[output:956c781d]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:28892c8b]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.06 rmse = 66\/0 muV ae = 29.2%)\nAllocating memory...\n","truncated":false}}
%---
%[output:10c20afb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:43df5a9d]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.06 rmse = 47\/0 muV ae = 104.8%)\n","truncated":false}}
%---
%[output:0a3e506b]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:4c23c1a8]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.64 rmse = 66\/0 muV ae = 19.6%)\nAllocating memory...\n","truncated":false}}
%---
%[output:611a36b2]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:4032f9d3]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.64 rmse = 45\/0 muV ae = 69.2%)\n","truncated":false}}
%---
%[output:4efbe68f]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:8fcde6f3]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=-0.04 rmse = 66\/0 muV ae = 16.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:480a9f12]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:9cc25b1c]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=-0.02 rmse = 47\/0 muV ae = 158.0%)\n","truncated":false}}
%---
%[output:734af950]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:00562850]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=0.99 rmse = 0\/0 muV ae = 0.5%)\nFBAR (r=-0.00 rmse = 67\/0 muV ae = 28.0%)\nAllocating memory...\n","truncated":false}}
%---
%[output:128cbb43]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:8736beca]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.00 rmse = 47\/0 muV ae = 102.0%)\n","truncated":false}}
%---
%[output:0f211e6e]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:77c78d20]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.06 rmse = 66\/0 muV ae = 14.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:99b43421]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:19e4bff7]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.04 rmse = 47\/0 muV ae = 67.3%)\n","truncated":false}}
%---
%[output:93482ee6]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:973ec693]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 1\/0 muV ae = 11.9%)\nFBAR (r=0.98 rmse = 66\/0 muV ae = 17.1%)\nAllocating memory...\n","truncated":false}}
%---
%[output:5923081a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:44dc52b1]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.95 rmse = 174\/0 muV ae = 203.2%)\n","truncated":false}}
%---
%[output:338c7722]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:6afd1434]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.4%)\nFBAR (r=0.10 rmse = 66\/0 muV ae = 19.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:2e1b636a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:48a354c4]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.09 rmse = 47\/0 muV ae = 63.5%)\n","truncated":false}}
%---
%[output:565458d5]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:5d8e6577]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=0.99 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.05 rmse = 66\/0 muV ae = 24.3%)\nAllocating memory...\n","truncated":false}}
%---
%[output:572d0ab0]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:8d42e345]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.04 rmse = 47\/0 muV ae = 50.2%)\n","truncated":false}}
%---
%[output:01fc80a2]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:4869747d]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.11 rmse = 66\/0 muV ae = 19.4%)\nAllocating memory...\n","truncated":false}}
%---
%[output:08b15ee4]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:237161ab]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.11 rmse = 47\/0 muV ae = 57.1%)\n","truncated":false}}
%---
%[output:5a393e9a]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:918abdb6]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=-0.03 rmse = 66\/0 muV ae = 33.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:90728558]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:5125bbce]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=-0.02 rmse = 47\/0 muV ae = 234.4%)\n","truncated":false}}
%---
%[output:72713035]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:2fe8c8a1]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=0.99 rmse = 0\/0 muV ae = 0.4%)\nFBAR (r=0.05 rmse = 66\/0 muV ae = 21.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:5d2383e2]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:1c7d796a]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.05 rmse = 47\/0 muV ae = 55.0%)\n","truncated":false}}
%---
%[output:8927da9d]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:4d18e78a]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.05 rmse = 66\/0 muV ae = 25.0%)\nAllocating memory...\n","truncated":false}}
%---
%[output:76589a7a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:5180fd7c]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.05 rmse = 47\/0 muV ae = 153.7%)\n","truncated":false}}
%---
%[output:6ead7dc4]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:2e974166]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.10 rmse = 66\/0 muV ae = 26.4%)\nAllocating memory...\n","truncated":false}}
%---
%[output:468d7a24]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:9da770da]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.08 rmse = 46\/0 muV ae = 146.0%)\n","truncated":false}}
%---
%[output:396b2734]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:83190844]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.02 rmse = 66\/0 muV ae = 34.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:75e301c5]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:4af4b559]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.02 rmse = 47\/0 muV ae = 137.7%)\n","truncated":false}}
%---
%[output:46087d98]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of RAW data.\n","truncated":false}}
%---
%[output:3105f2b4]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.10 rmse = 66\/0 muV ae = 8.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:23f9252a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:72f65666]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.10 rmse = 47\/0 muV ae = 88.0%)\n","truncated":false}}
%---
%[output:70ab0732]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQnYVVX1\/5cDYuKAM0YDZkI4xKAipuaQWuCAOeeAUc6mUA4gDgiKaIWAQ2UGqSkqCmKEYA5pqJkTmRaCpJai\/tRwwpSY\/v\/P\/rnf3+Fy793r3Peec889Z+3n4VF499nDd6313d+99j7nXW3FihUrxIohYAgYAoaAIWAIGAJNhMBqJmCayFo2VEPAEDAEDAFDwBBwCJiAMUcwBAwBQ8AQMAQMgaZDwARM05nMBmwIGAKGgCFgCBgCJmDMBwwBQ8AQMAQMAUOg6RAwAdN0JrMBGwKGgCFgCBgChoAJGPMBQ8AQMAQMAUPAEGg6BEzANJ3JbMCGgCFgCBgChoAhYALGfMAQMAQMAUPAEDAEmg4BEzBNZzIbsCFgCBgChoAhYAiYgDEfMAQMAUPAEDAEDIGmQ8AETNOZzAZsCBgChoAhYAgYAiZgzAcMAUPAEDAEDAFDoOkQMAFTxWTdunWTDz\/80NVYc801ZbPNNpMDDjhAzjnnHFljjTXcv5fW2WSTTWSvvfaSCy64QD7zmc+0tP7nP\/9ZvvOd77i\/77777nLjjTe2\/Ozvf\/+7a9eXtdZaS77whS\/ISSedJIcddthKI7z55pvl17\/+tTzwwANVnW3p0qXSuXPnVepccsklcswxx6gddcaMGTJ69Gh59dVX5bOf\/aycfvrpq4wp2tj7778vw4cPd+NbtmyZfP3rX3d\/33TTTV21N998Uy699FIBD7AFr29\/+9vywx\/+UFZffXXx81uwYIGAw1ZbbSWDBg2SPfbYQ\/i9o+A6bdo0WXvtteVHP\/qRHHXUUa7d2bNny\/HHHy\/333+\/s1OlMnfuXLnsssvkL3\/5iyxZskS22GIL6dmzp\/zkJz9xj4A5bVDatm0rG2+8sXTv3l1OPvlk2X777dW4UXHOnDlu7oxt\/fXXd\/M899xznS+VK9XmTn3scO211670KPg98cQTdcEm1uSscgsCUZ\/BR\/H1Xr16ybBhw1r8PlqHBzfaaCPp0aOHDB06VLbccsuWtj7++GPZaaed5D\/\/+Y+st956zrb4oS8aTqLuiy++KPvtt5973sdeJZM98sgjMm7cOPnrX\/\/qYgKeGjlyZEt15sTPb731ViG+u3btKhdddJEbf9wyYcIEF\/8U2jzwwAPLNvHGG2\/Ieeed52IHngCjKOfddtttcuWVV8onn3zi2qDN1VZbTT744APZZ599XIzz30qFeTDHhx56yD1DHDGviy++WDp27CjRccJD2KJLly5y9NFHS9++feNO29WvNHf+fdKkSTJv3rxVcHnnnXccNz711FPyuc99Tn784x+34P7zn\/9c7r33Xrnrrrvc3ItYTMBUsbonC\/4LiUAGlKgI8HW+8Y1vOCdn8Vu0aJGceeaZbuH1hYBngaIgfv70pz+5oKFEBQzih+B94YUX3M9mzpy5khA58cQTnbi58MILq\/qrFzAbbrihbLfddi11WeT33ntvla\/7cRHAe+65pzz22GNubszja1\/7Wtk2WOjvu+8++cpXvuJEBkIBQr799ttdff9zCH7bbbeVKVOmOFJkYe7Tp48j\/X\/84x9OhNA\/QU3\/YP\/ss886kUL\/Tz75pFxzzTWO4NZZZx0nDvbff38nQCqV5cuXO0H1+uuvO2G0+eaby7\/+9S+BJOiL4hcahE379u3l5ZdfdiTZpk0b+cUvfuHEqaYsXrzY9fX22287wfrSSy8JogwBeNZZZ5VtotrcEUBewOBz\/J2ywQYbyFVXXSWzZs1qFTaaOVmd8gh4n+nUqZN8\/vOfl2eeeUY++ugj549XX331Sn6FzxOP2As\/JA5YjH2555575Ac\/+EHL3\/E5hIgvGk6i7vjx411sTZ8+PWi2m266yS2gCBVEd6mAYcME5yGEvvrVrzohAdf98Y9\/dP6nLcQSMc7Ghj\/VBAzjIFbYNME7UQGDyGMc\/Bz8jj32WLepY5PDOIk1\/l6tsPmZOnWq42C4ijhlc4MYAGMvNtZdd13HFWzgFi5c6JqsFsOV+qw2d8aCL8BFbPCiuPz0pz91fMc6wOYHMccYqYdfsBGuRUhqbZb1eiZgqljIk8XEiROld+\/ecvbZZztSQIX7XURpHZT\/r371K5dRYWGhICZ4ngAgWBAnI0aMcIFHiQoYgg8iYfFjwYMAIULKf\/\/7X5ctQHmzKOLwjAPCJHtBACNsWNy8gCnN9kSn+\/3vf1\/+8Ic\/VCQSH+TsAM444wyB6NihMLYbbrhhFeQIUsRRu3btnMBg0d91110dOYAbmQx2L8yfthkbO1DIymMcbZSdEc9QEIavvPKKnHDCCUJWiB0J2Rh2jZD09ddf74KcPisV8KRPSIsMkN+1IJggKYpfjLwAhTAQouAEmULaYB3CDsHGDhK7MzeE2Le+9S2BEBk7oqxaKZ37l770pRYBM3ny5FVI68EHH2wVNlknqiyPr9Rn7rzzTrfYYDOfzSutw7\/zbz6D5ucX3QAQJ8QLQt0XDSdRt3\/\/\/m6Rh7PYdJCt+P3vf++yDbTBhmrrrbdeCVZ4hUxkqYBhs8KCSQZm5513dgILocUcTznlFBd7o0aNWuW5aONsHnw2GX549NFHqwoY\/ywblEMPPXQlAcN8yIbCRbvttpvsu+++8stf\/tLVIRsDH0SzWuV8h+cQDdgBO1HeeustlzVHnHkB4wUm40dYeEGKiKOPes6dzRnCNipgWE8QLIg4OJP1gfnByWRyfeY4y\/GR5NhMwCgFDLsmFk8yAQQOBEEpzcDg2JAEpOOFx8MPPywDBgwQdmgEP6QSzUqUChgWWp6lnWgGhlQvGRiyGgiUb37zmy4IjzjiCCdu2FEgYth9eAHDQokg4iiEQOf4i0WUElqE2S2xK0GQIUw8mbATQwCUFoQFuxOIk7FQfFASiBz3sLAjWhgThf+CJ\/V8YZ4QHP2BN9ktiII5nXbaaU5EUMDSZ5TYqTD3aoVd8Y477ihkRxATkDHkHCXy0oWG9hBmjIECeZBqDmFHNuU3v\/mN8xnmS0G8YidPfuXGWmnu1PUZGAiWVD8pbS8EW4tNkiST97ajPoNPYqfrrrvOZRv8kZ+v4zMw2Jk4P+6449wxIwWxDC\/gnwh+Fm7iF8HrYzYqYCpxEhkKduX4H+0RM\/AImwn6Z\/OBDyHKfbv0X07AvPvuu7LDDju48f3tb39zC7yvh1hgsdUs4mSSqIvwId5LF+pKPlJOwFD3Zz\/7WcsGkQ0VOBOTzG\/w4MFBl4OL4BbwJLbhAjaH\/mpAqYDxXIXwIUNOfMM99Zx7OQGDj7B2vPbaa+4UYOzYsc5m3\/ve95z4Ch0PBoFo8gomYKoYMHre7KuxO4Fw\/D2GcnVwKgjMZw985obFjAWeRZQFBxLjvLX0Dgx9QVz0c+SRR7aMEBFAtoD0MAspWRF2CAQz5aCDDnKChnNdsgUszvRFULJochQCKXrVTv\/srAheUt+lhbY5XvHZEcQMpEx7nLGXlltuucVlgKJpcU+eLLSQO2KAjAbnuWRDIHuEBYLH39mJ3vXgeIhMRvTeDgQN\/mRbyGSR3qV98CLbw6610jENIovd57\/\/\/e+W4UPs7OAg53IChnT3l7\/8ZZexIQsFiYWwwzbYKHqUyKLCcRkijoWlXKk29zFjxrh+EcLYAr\/BT5iT33G2Bpsm57KGDb\/0fgsDIe5YVP0CU64OcURMH3744W7sPnPjNwBsYshOkj05+OCDXR0NJ7Gw4f9PP\/20IECIb+KIo10WQQTE7373O5c1ifJLOQED37DxYaNB7FI4tiDW4BeON8gukgnFB73YiRqDn8NNCHk2fuUW6rgChvpwKEKeuCW7BPcg\/C6\/\/HI3HrAiQ+2PW0vHRBYDfH354he\/6DZrZGPLCRjqgRfZZX+MVM+5V8OFzSyiEzuAJUfmHLNjMwrj8RvmhgVCAzo2AVMFdE8WLLSo7vnz57vF4+6773bOFCUUFnkWJRY4zmH9mS27Kf6dtKdfuLz6Z2EmoxIVMOxqSBeywJamjznz5NgJEmBHw4WucsWfB5P25LiD4rMjBDMZHE3xGRgEE3c\/tBkYcCPtSSnNwLDbgQjZVXKR0S\/YpM6jOyeIiXGy++BCoz+CKiWhQw45xO3qSGUjwugPjBAklS7xgctzzz3nsjzgiG18hkibgQnh5zMwtDdkyBBXnYwJ86qWgaFepblH7cn\/c0zJMQPEDU71wCY0L\/v5qgh4nyHDxr2pxx9\/3N3\/4siFWKCU+hX3ufBdCtkINhzENQLVH9my2+YYmtgjBqN8U42T8D02HggoYsj3Uzry0pirNQMT8gmOXRDfZEjhI+bOcTqZTOYcFVGlbVXKwETrwbHEOqKNthFqZEbYMLKo+wxouXFyLA3fIlg4nunXr58ba6UMDJsdNn0+A1PPuWuEHZtE+J11Br\/AV\/zxFvcq4dQiFRMwVawdTdciQlC+KPZo4JeeSbPrYUfF7gp1TuoWdVyubLPNNm4nVHqERBaFhZ7A9GLE398gu8IlXt7EGThwoNvx8F9fOKJgB8fuhFSjTxF7AcOlXsZIgWjpix0amaDSQnAg1nz2xN+B8fdqWGgZB9kQsCl3B2aXXXZxYox7G4zLH9dA1JC2vzPkMSVD4S8Gkvng8jBzj+5C\/TgRKmQ1GB8ZEjIrkADn47QXvQzJMwgV+oXsfAbNZ4gQGSwy2jswIey4mAlxRnepHFv5N0sg8ih2jC80dwjWn9ezE4OcETDRI81asSkS6dV7rqU+Q3aCLEVUyJfWIXb4OVlRNjZsjLgv5Y9Wo2PEVzmyJXY1nMSRCpyDMCAjSXzDAwia6NsqZIcQ1b5UugPDuLgf4jdgtA2fcBx96qmnumwsGwJ802edo+PnOJ34LVf8ESt3uN577z0X79FFWCNgEEiIQC4ikyGHP9iYsDmEQ+Gt0sL4ESM+O+OzSv4IvlTAwEWIydI7MPWYux9bSMDADxzlcxRH5pWNMFlexgb\/3nHHHWUzYPX29yy1ZwJGKWAIYsQICx5kwO6dxSh6B4YLef5sm0Ama+CDnbRk9GIZuwTIhN04aX\/\/GjWLFMWToD+OIcPDLsy\/SshizIKIAIGoEAP\/\/Oc\/3c9Jcf72t791CxvkRYoVgiCTET1zD93jeP75592lOHaTCCpIgnN6f4ziL5oSTP6tqehZP88hlhBZBBeF3SC7Qv6NNDvn9IzLv9kEnmAHuSIW+UM74BUVWRAQ6WHS5cyP8SGQIC0W9ugdJG9iT+aQFkIKkceOjUWDjBH9+vGTnoVIEWWQIMdVZHX8PZsQdtiUhQTxBulAdBx1+dRzOexCc2fHxQVIjtpI7bNocDQAFtEjwFqwyRIpNdtYSsUJiz22ZzPhj1+jcYGfsmnB9\/BtuAQxS7xyVy36VgkX9MlW4Ou8PFC6YSrlJDIvxAJZhQ4dOjgo\/cVgjj7ZiOGT3CPj+JV\/4y4Iiz8xzLjgKcaA35IF9os57bE54LgG7uNuH6JKcw8katNyC7W\/3O8vqBOrV1xxhRsr\/cC5ZKA55uF41hc2N\/Ag4ooNIZsGNhfwCQs8MePvGEXHAEY8Cy8jeOApRKXPivs5E28IM7iVmKVE30Kqx9zxEewMvzJvbEQ8wzFkqXzBP8j+IM7gfTaS2BABQ9YfP6r2CYlmiyvNeE3AxBAwLHScByMy\/O4jeibtL8uyMJ9\/\/vkuuPl\/iIxMRvQ7Ij6IuSNBYJYKGIKL4GORJajZHbGAk7r0hcWV816cH2HBzyEkHJ3jLlLQBCY7G14ZRqWTlfDflQgtwvTD8Yz\/DgwE5nd2\/KzcIsyZO\/1zQZCxE2QQiCdTLqNx9g5p8jzZJHZhPo0MOfEz2oFYICXGHL0zgqCAgDhy8t+RYIcKebF4gCXfePDHZx4vnmP3AlFASBAWZM3i4u8YRO8qIMwQpQgbjqiir6NrsEN8gQULFQRMeppMD2KoHHahuUOW3BXi7TPmhj+RNo8ueLVioyELq1MegXLHjsQ\/R0g+W1l6BwZ\/QERjP7J0ZG2JVZ8J9D35rIjfyJQKmFJOIra5S0Pc+oKvwQUIDwQOwhxeIn7gDJ8tLJ2dz4r6Iwr\/HRhikjcA\/X2XeizipQImenE+Oq5oVot\/hxvgCf\/dGhZ47tgRc8QtWZNyF12ZM5s8eBJ8qAMH8CzZ8+g3W4hXNj1kq8jwRO+a1GPuPtNdir\/PvvPv3Hnj\/iKClfuDFDaSPiPEuON83ysvsVwYAYPDksHAWdnlIATITlBYLHEAnJkzbM5Ay30ErlFGZ6HlhjxZhdCbNo0ao\/VrCBgCjUeAjRECQ\/MmTuNHayMwBFqHQGEEDEco7FhR2uzAESwctbCTJUOAymYHz5ki6pusR1YKaUUub7Ejin7dNyvjs3EYAoZANhAgY8Mdr9JvvGRjdDYKQ6C+CBRGwERh40iH4wCOXkjlkp7lDJiLdKRLyXaQqvPHHvWF3FozBAwBQ8AQMAQMgdYiUEgBQzaDS3McK\/kPvnHG6L\/iynksZ6pF\/kRzax3LnjcEDAFDwBAwBJJEoHAChtv5XLDjljoX2PwXTaMChu+fcNG03CuBSRrD2jYEDAFDwBAwBAwBHQKFEjDceeHWP5\/S9t\/TAKbQERKvL\/JqnhVDoFkR8L+TqVnH38hxW\/w3En3rux4I5DX+CyNgeK2OV4B57S36bj3OwRtJvH7oL\/H63+TqHQex47\/PUs6ZeDU3+jG5ejhcEm00yziZe7OMtVnGGfLhJPwtL22GsOPn\/a6YUXG6dw\/uU5U\/0sKpWXzV4r\/+HhHy4fr3mE6LhREwfLKaz9dHfy8Gd2C4uMtr1AgQ\/xo13z2Jipy8Gj8dF7NesoCA+XDtVghh1ywCpnYE7MlmRyDkw806v8IImNYYKK\/Gbw0m9mxzIWA+XLu9QtiZgKkdW3syHQRCPpzOKOrfiwkYBaZ5Nb5i6lYlJwiYD9duyBB2JmBqx9aeTAeBkA+nM4r692ICRoFpXo2vmLpVyQkC5sO1GzKEnQmY2rG1J9NBIOTD6Yyi\/r2YgFFgmlfjK6ZuVXKCQNF9uDW\/SiSEnQmYnARJjqcR8uFmnboJGIXl8mp8xdStSk4QKLoPt+ZXiYSwMwGTkyDJ8TRCPtysUzcBo7BcXo2vmLpVyQkC5sP\/Z8i4v0okhJ0JmJwESY6nEfLhZp26CRiF5ULGD33oKq8fEVJAZ1UygkDIhzMyzFSGEfdXiYSwMwGTitmsk1YgEPLhVjTd0EdNwCjgDxnfCEwBolVpKAIhH27o4FLsvJZfJQJ2fCeq0scqLf5TNKB1FRsBPmDIn2ofY43daEYeMAGjMESI\/I3AFCBalYYiEPLhhg4upc5r\/VUiIews\/lMyoHVTMwIhH6654QY\/aAJGYYCQ8Y3AFCBalYYiEPLhhg4uhc6T\/FUiFv8pGNC6aBUCeY1\/EzAKtwgZ3whMAaJVaSgCIR9u6OBS6DzJXyVi8Z+CAa2LViGQ1\/g3AaNwi5DxjcAUIFqVhiIQ8uGGDi7jnYews\/jPuAFteBLy4WaFyASMwnIh4xuBKUC0Kg1FIOTDDR1cxjsPYWfxn3ED2vBMwBTZB4zAimz9fMw95MP5mGUyswhhZwImGdyt1fohEPLh+vWUbkuWgVHgHTK+EZgCRKvSUARCPtzQwWW88xB2Fv8ZN6ANzzIwRfYBI7AiWz8fcw\/5cD5mmcwsQtiZgEkGd2u1fgiEfLh+PaXbkmVgFHiHjG8EpgDRqjQUgZAPN3RwGe88hJ3Ff8YNaMOzDEyRfcAIrMjWz8fcQz6cj1kmM4sQdiZgksHdWq0fAiEfrl9P6bZkGRgF3iHjG4EpQLQqDUUg5MMNHVzGOw9hZ\/GfcQPa8CwDU2QfMAIrsvXzMfeQD+djlsnMIoSdCZhkcLdW64dAyIfr11O6LVkGRoF3yPhGYAoQrUpDEQj5cEMHl\/HOQ9hZ\/GfcgDY8y8AU2QeMwIps\/XzMPeTD+ZhlMrMIYWcCJhncrdX6IRDy4fr1lG5LloFR4B0yvhGYAkSr0lAEQj7c0MFlvPMQdhb\/GTegDc8yMEX2ASOwIls\/H3MP+XA+ZpnMLELYmYBJBndrtX4IhHy4fj2l21JhMjCzZ8+WESNGyJw5c2TAgAEyePDgFqTfffddGTRokMyfP1\/at28vY8aMkc6dO7f8PGR8I7B0ndZ6i49AyIfjt1icJ0LYWfwXxxeadaYhH27WeRVGwCxYsEAWLlwo9957ryxbtmwlATN8+HBp06aNDB06VKZPny4TJkyQyZMnm4BpVq+2ca+CQF4JLA1Th7AzAZOGFayP1iAQ8uHWtN3IZwsjYDzIY8eOlcWLF68kYHr16iWTJk2STp06yfLly6Vnz54yc+ZM6dChg3ssZHwjsEa6sPWtQSDkw5o2ilonhJ3Ff1E9o3nmHfLh5pnJyiMtvIBZunSpOy6aO3euy8JQ+vbtKyNHjpQePXqYgGlWz7Zxr4RAXgksDTOHsDMBk4YVrI\/WIBDy4da03chnCy9glixZIl26dFlJwPTp00dGjRol3bt3NwHTSO+0vuuGQF4JrG4AVWkohJ0JmDSsYH20BoGQD7em7UY+W3gBA\/iaI6SBAwcKf8oVI7BGurD1HUJg3Lhxwp+XXnopVDW3P7dL\/Lk1rU1MgYAJGAVIzVCl3B2YYcOGSdu2bVsu8Y4fP16mTJnSMp2Q8U3ANIPliz3GkA\/nHR27xJ93C9v8qiGQ1\/gvTAZm3rx50r9\/f\/noo4+cndu1ayfXXXeddOvWTXiNmuyKf4169OjR0rVrVxMwxgm5QSCvBBbXQHaJPy5iVj8PCOQ1\/gsjYFrjhCHjazIwof579+4tEydODFWznxsCNSEQ8uGaGm3Ch0oFjF3ib0Ij2pBjI5DX+DcBo3CFkPE1AqbfFTOq9nT34D6FvqOgMINVaQUCIR9uRdNN9WipgLFL\/E1lPhtsjQjkNf5NwCgcImR8EzAKEK1KQxEI+XBDB5di57UeIdkl\/hSNZF3VFYE8X+I3AaNwlRD5m4BRgGhVGopAyIcbOrgUO7dL\/CmCbV1lBoG8xr8JGIWLhYxvAkYBolVpKAIhH27o4FLo3C7xpwCydZFZBPIa\/yZgFC4XMr4JGAWIVqWhCIR8uKGDy3jnIew08V\/kb\/Bk3LyFGF7Ih5sVBBMwCsuFjK8hMLvEqwDaqiSGQMiHE+s4Bw2HsNPEvwmYHDhCE08h5MPNOjUTMArLhYyvITATMAqgrUpiCIR8OLGOc9BwCDtN\/JuAyYEjNPEUQj7crFMzAaOwXMj4GgIzAaMA2qokhkDIhxPrOAcNh7DTxL8JmBw4QhNPIeTDzTo1EzAKy4WMryEwEzAKoK1KYgiEfDixjnPQcAg7TfybgMmBIzTxFEI+3KxTMwGjsFzI+BoCMwGjANqqJIZAyIcT6zgHDYew08S\/CZgcOEITTyHkw806NRMwCsuFjK8hMBMwCqCtSmIIhHw4sY5z0HAIO038m4DJgSM08RRCPtysUzMBo7BcyPgaAjMBowDaqiSGQMiHE+s4Bw2HsNPEvwmYHDhCE08h5MPNOjUTMArLhYyvITATMAqgrUpiCIR8OLGOc9BwCDtN\/JuAyYEjNPEUQj7crFMzAaOwXMj4GgIzAaMA2qokhkDIhxPrOAcNh7DTxL8JmBw4QhNPIeTDzTo1EzAKy4WMryEwEzAKoK1KYgiEfDixjnPQcAg7TfybgMmBIzTxFEI+3KxTMwGjsFzI+BoCMwGjANqqJIZAyIcT6zgHDYew08S\/CZgcOEITTyHkw806NRMwCsuFjK8hMBMwCqCtSmIIhHw4sY5z0HAIO038m4DJgSM08RRCPtysUzMBo7BcyPgaAjMBowDaqiSGQMiHE+s4Bw2HsNPEvwmYHDhCE08h5MPNOjUTMArLhYyvITATMAqgrUpiCIR8OLGOc9BwCDtN\/JuAyYEjNPEUQj7crFMzAaOwXMj4GgIzAaMA2qokhkDIhxPrOAcNh7DTxL8JmBw4QhNPIeTDzTo1EzAKy4WMryEwEzAKoK1KYgiEfDixjnPQcAg7TfybgMmBIzTxFEI+3KxTMwGjsFzI+BoCMwGjANqqJIZAyIcT6zgHDYew08S\/CZgcOEITTyHkw806NRMwCsuFjK8hMBMwCqCtSmIIhHw4sY5z0HAIO038m4DJgSM08RRCPtysUzMBIyLvvvuuDBo0SObPny\/t27eXMWPGSOfOnVtsGjK+hsBMwDRriORj3CEfzscsa5tFGvFvAqY229hT9UEgr\/FvAkZEhg8fLm3atJGhQ4fK9OnTZcKECTJ58mQTMPWJHWulDggcffTR8vjjj1dsqXfv3jJx4sSKP88rgdUB2lTi3wRMPSxlbdSKQF7j3wSMiPTq1UsmTZoknTp1kuXLl0vPnj1l5syZ0qFDB+cvIeNbBqbWsLLntAhofKzaIhnyYe048lgvjfgP4RYSoKHn7efxEQhtCmgxK3YJjTU0zrzGf+EFzNKlS91x0dy5c10WhtK3b18ZOXKk9OjRI1cCZty4cTJw4MD4kd6AJ5plrGmN0wRMMk5o8Z8Mrq1tNY24CsUUc7h7cB+ptjFIY5x+E13tGkJonCZgWuuRGX1+yZIl0qVLl5UETJ8+fWTUqFHSvXt3N+qQ+s3o1GxYhkALAqEdWlGhsvgvquWLNe+8xn\/hMzC4cSiFXCxXt9kaAsVCwOK\/WPa22eYHARMwIjJs2DBp27ZtyyXe8ePHy5QpU\/JjZZuJIWAIVETA4t+cwxBoTgRMwHz6GjV3Q\/xr1KNHj5auXbs2p0Vt1IaAIRALAV6jtviPBZlVNgQygYAJmFaYgSzNlVdeKStWrJB+\/frJueee24rW6vMoY7nkkkvcW1RrrLGGDBkyRPbff\/9VGn\/99ddljz32kI033tj9jAzUww8\/XJ9B1NjKBx98IKeddpo899xz7s2vu+66q8aW6vvY7NmzZcSIETJnzhwZMGCADB48uGwHWcS0vkhYa1EEshj\/fnzEPRyw9tpru38irvr37+\/+\/\/nnn5cEikuCAAAgAElEQVRzzjlH3n\/\/fdl2223dd6\/WXXfd1Ix7wQUXyIMPPihvvvmmiym4x5dKmIa+1ZPU4CuNtVqsN2qsSWGQ5XZNwNRonbffftsJg6lTpzoRcMghhzixsPvuu9fYYn0eg7Q4AuObIG+88YYceuihjizWW2+9lTogAI877jh54IEH6tNxHVpZvHixIBYWLVok1157bWYEzIIFC2ThwoVy7733yrJly6oKmKxhWgezWBNlEMhq\/EcFzK677ioHHnjgKqPnLcuzzz5b9t57b7noootk\/fXXd39PqzzxxBOy5ZZbys4777ySgKmGaehbXUmNvdJYq\/Fno8aaFAZZbtcETI3WQSA89dRTLgND4eN38+bNk8svv7zGFuvz2JlnnikQ15FHHukaPOmkk9xr4QcffHDmBYwfIB9su+KKKzIjYPy4xo4dK4isahkYEzD18eOst5LV+A8JmJdfflmOPfZYefTRR11VMiCnn3662+SkXciyRjMw1TBt9EXr0rFWEzCNHmvadmxkfyZgakSfezIc1\/idC5mPW2+9VW688cYaW6zPY0cccYScccYZLZkgdgObbbaZnHrqqasImL322ks233xzWWeddVw25phjjqnPIFrZSjMLmKxi2kqT2OMlCGQ1\/qMC5qGHHpI111xTtttuO7nwwgulY8eO8sgjj7gjI\/+lcY6RyIS88MILqdu4VBRUwpSMcuhbXUkPvpyAKRfrmu8KJT3WIrVvAqZGaxNslLPOOsv9d8aMGXL77bfLDTfcUGOL9Xns8MMPdxcSd9ttN9cgb1hsscUWcsopp6zUAZkEjkX4Gbsyzscvu+yyhh+BMchmFTBZxrQ+3mWteASyGv9+fBx7snGhcBw7a9YsJ1oQMHx87Y477nA\/477GLrvskhkBU45Tr7\/++uC3upL2zFIBUynW+d5K6LtiSY+1SO2bgKnR2llNIWuPkEqnzZeHOQsne9Po0qwCJsuYNtqmees\/q\/FfDmfulPHrUfja+CuvvOIyrY899pirakdIOs8sFTDVYt2OkHSY1qOWCZgaUfQXzu6++27ZaKON3GVZ7kY0+hIvmSDu4\/hLvFwu5nwbcfLss8+6y7wEIzu0DTbYwL19wNsAkNp5550n++yzT42I1O+xZhIwzYJp\/axjLYFAVuPfW+fFF1+Urbfe2h1zwwe8bDBt2jT3Y740zltI\/hIvHNCINyhLRUE1TBv9rZ7SsVbjz0aPtUgRagKmFda+8847hYud\/ALIgw46yL2F1OgCYfHKL2\/MrL766k6U+Neoua+z1VZbufsw3Nm59NJLhU+p8xojvy6h9JipEXPZd9993dEWu8YNN9xQTj75ZPfqciMLl7M5Yvvoo4\/cMNq1ayfXXXeddOvWzd2ByjqmjcQuz31nMf493lzUJeMCB3CkcfHFF7uNCwXRzWaLzxZss8027kUENjhpFXiS+zlvvfWWO+baaaed5Oqrr3bdV8K0Ud\/qqTTWavzZqLGmZb8s9WMCJkvWsLEYAoaAIWAIGAKGgAoBEzAqmKySIWAIGAKGgCFgCGQJARMwWbKGjcUQMAQMAUPAEDAEVAiYgFHBZJUMAUPAEDAEDAFDIEsImIDJkjVsLIaAIWAIGAKGgCGgQsAEjAomq2QIGAKGgCFgCBgCWULABEyWrGFjMQQMAUPAEDAEDAEVAiZgVDBZJUPAEDAEDAFDwBDIEgImYLJkDRuLIWAIGAKGgCFgCKgQMAGjgskqGQKGgCFgCBgChkCWEDABkyVr2FgMAUPAEDAEDAFDQIWACRgVTFbJEDAEDAFDwBAwBLKEgAmYLFnDxmIIGAKGgCFgCBgCKgRMwKhgskqGgCFgCBgChoAhkCUETMBkyRo2FkPAEDAEDAFDwBBQIWACRgWTVTIEDAFDwBAwBAyBLCFgAiaGNSZMmCCXXnrpKk+cd955cuKJJ0rpz9dbbz350pe+JIMGDZI99thjpeeOPvpoefzxx92\/TZ48WXr06NHy85NOOknuv\/9+9\/cVK1bIpptuKr169ZJhw4a5\/\/eFn\/HvF110kRx44IFVZ\/LII49I\/\/79V6lz3333yVZbbaVCgf7GjRsnt956q7z\/\/vvStWtX13d07NGGli5dKp07d16l7UsuuUSOOeYY9+9nn322w+H11193f\/\/zn\/\/cMsfZs2fL4MGD5V\/\/+pfsuOOOMmbMmJafHX\/88dKlSxcZOnRoxbEz3p\/\/\/Ody2223yVtvvSXYg7n+4Ac\/kN12203+\/ve\/ywEHHOCeX2ONNeQzn\/mMbLnllvKtb31Lvv\/978taa62lwsXbKQ420YZPPvlkwQ7R+Yeway026olZxSACUT+KVt5rr71k\/PjxK\/kZP2\/btq189rOfFTgAP4uWq6++2vk55ZxzzpFTTz215cel\/LLuuuu6GBwyZMgqMYiPb7DBBjJy5Mjg+KvFIA\/PmTNHhg8fLvjc+uuvL9\/+9rfl3HPPlTXXXDPYNhXw5R\/\/+Mcyffp0eeedd1wbX\/va1+Tiiy+WDTfcsGob0TkTX\/Acbfzwhz+Up556Sj73uc+5tj0HEe\/33nuv3HXXXbLaaqtVbPv3v\/+947KXXnrJxfnnP\/95OeSQQ+R73\/uee6Zbt27y4Ycfuv9v166dbLbZZrL77rvLaaed5v4\/blm0aJHjFXgOjnnggQdcEyFehruuvPJK+eSTT9zcWX+Y1wcffCD77LOPXHbZZe6\/RS0mYGJY3gcTxMHi6QvCAOfyP99kk01k3333dcT17LPPOsJ65pln3AJJefvtt6V3795OnFC++93vOiHgixcwnTp1coHFsx999JHsv\/\/+AsH58vzzz0u\/fv1cIIeIwAcKwUPQ+zJq1ChHppry61\/\/WhAfiKivfvWrLggRBX\/84x8dWZYWvwgztu22267lx4iPvffe2\/39G9\/4hhN5Dz74oMMjKmAOPfRQWWeddWT06NHyzW9+U4466ignaBB3iEb6hwwrldtvv93VW3vttaVnz54CicydO1d+9KMfOcEZXXh22GEHeeONN1qE1M477yy\/+c1v1CQdFxs\/5ilTpjgRx9whJj\/\/EHatxUZjb6ujQ6DUj\/xT+BTiIvrz73znOy7+EazYe+LEiY4LfIE3\/vGPf7i\/fuUrX5F77rmn5WeeX\/B5Fux58+Y5n914443lySefbKm3bNkyoe+f\/vSnqsWtWgwuXrxYvv71r7sxs4Cz4C9YsEBOP\/10Oeuss1QA+dj44he\/6HjxoYcekvnz58txxx3nhFGl8vLLL0ufPn2E+fDHCxjmdfPNN8vMmTOdkEJoIFjefPNN2W+\/\/eTGG2+suKmiL+aAmCDGunfvLm3atHGYI1oQnBQvYNiA0fcrr7zi\/gv3EbMdO3ZUzd1Xwg\/YqNJGOQFTjpfhTXgWrNmoHnvssQKWbIbhYebB34tcTMDEsL4nEJwJZVxaSn\/OgokDUlhscVKKD2gI6oUXXnAExMK1+uqru597AXPmmWe67M2dd97pApWF3mdmqHfNNdfIH\/7wBxcYLIAsuPyBYBA\/LNRenXsBU7qr83NAaLGzigZX6fzYNUESZGBY4NnlQbCM7ZRTTlkFD78IQ3yQSrUCFv\/9739XEjCIM4jiqquuciQKmdEXJDVw4EC3Y6pWEAaQTTTjgxB87733XLvRhQUyoDz88MNywgknOKK54oor5PDDD3citN7Y0BdYIsxY1G644YaV5h\/CrrXYxHB7qxpAoJwfRR8p93M2PH\/729\/cjppMDIVMB5sU+ABxQ6aBbMLWW2\/tfl7KL8Q5sUWJCv+nn37atcl\/2Ww99thjTsy8+OKLbqPDRoB49Xzjx1ouBv0mAJGF2EI0sfjTLhsnshfbbLONyxBExxCd\/4gRI5x\/E79klPh\/4tLzWzl4ly9fLocddpj7ERmQRx99tEXAkHVAsDAvOI7YJbtDVoas0E9+8pOqFvN8Gt0Q0t8\/\/\/nPFo72AsYLTMQUG1UwhwvYVFFCc6cOHE1mB7689tprywqYcrzM+rH99tu7TBUZY8TfL3\/5S\/c8\/sOc\/ZpS1CA1ARPD8p5ACKgvf\/nLLU8SUKRyK2VgIIZp06a5YwoKAcCiSGqQXQWBw44CgUCJChicnmC57rrrHAEQAL4cccQRbmGnDgKB3Yx37jvuuMNlEyAKdmtewDD2JUuWuAWcYxyfMg0t0u+++67b1VEgXrJJpGshC4KJeZQWvwhDcggsiJkgJFghwGgpR56ICQjq448\/duOF+CAuxCDzqpYipu2xY8c68UOGiTEy\/l122cURIqXSwgMR3n333U4o\/eIXvwgKmFqwof8BAwa4HTR9QZhRARfCrrXYxHB7qxpAIOpH2NEX4otFuFIGhhhgEfrCF77gHkEwE+fENWXSpEkrZTqiAoZNBAsjYpu44hjW8wt8wXEPnIIoOuigg1wWkowHmRqEDRkBeCYUgxxbsymiH39c62PVb8pCizjzRzSxIPuCeCA+\/ZhLISbu4BQ2SCzgs2bNahEwiAhi57XXXnPZbdoBS7iMDV70mL2c6cAKgQcmcDG8gBCMPlcqYGiHTSsYkAH7y1\/+4poOzZ2jHjiP48QjjzxSyJyWy8BU4uWf\/exnjsMocD38z7Hjtttu67LRRS8mYGJ4QKU7MOxSdtppp1XuwPimISREDjueV1991aUA2SmwgyFIycjg3BznUKJ3YHwb3AHBeX2QkTblWITdBMFGpoWdyK9+9SuXFv3d737nAv\/ggw92QgkBw99xfDIQHPsgANiZkclgEYaQIAJ2WKWFFCuBiBBhN0LxognhBVmWFhZhfsbYISraZ6dGEJfuksoJGNqjP+YKabDYIyroi50pf8ADkkVAlhbIA7Hk75fwc+bHQoEYrCRgWADAmnGziCSBDT5z4YUXOiHG8Vrp\/DXYtQabGG5vVQMIVLoDg\/hmc1Hp58Tt9ddfLxw5U9hls+kghskIwAMcuyBUKOX4p0OHDs5Xo\/fQEOuIFo5JL7jgApc54e4GCzxZP+59bbHFFi6rES3lYvCMM85wIiuaLaEv7sB53ps6dao7jkGU+GPyaLsIl\/PPP99lgBjDLbfc4jYFPsNZCi9ZHsaPWCDrwZFzVMD4+sQ3R9jEAfURI9xPYWNF4eiFMZUrxDg4c0TmC1kpMryUcgKGjA\/HOBTmAqeF5o4PIBpnzJjhniknYKrxMn2BLZtOsOXujueNyy+\/XP7617+6sZLJq3acntcgNgETw7Jxj5A4N0aYcH7qz29R1IiGXXfd1e1s\/G6AQMTRyVZ4AQOhtG\/f3tVht8Cuy+\/wIBUc\/4knnnBCxJNP6XSoT9YCQoymjP3OClHCri9Uas0yRPsliCGV6A7G91tJwETHBZmyU9lzzz1dKho8OLOGlMlwVSoLFy50GCJGEG6QHH+vJGA4nqI9n4FJAht2w\/SP2KT4RQq\/4OgLu2mx4\/lasQnNzX4eRiDuERIZRRZKYpjjQy7akhXhuJI7X\/w\/tmeDwgJL\/OIP0QwvfsOizs\/hAX9BH87hiBtxz\/0N7tfh8+WKz6RWi0HPE3ASWRsK9\/9YUKPH4tVQQrwQq9zfQ1D4bLDnptJn\/UVmNnpwFmKHGGaTwjzh1GhBELEJvOmmm1ymgwwq+MG5f\/rTn2SjjTYqOzzswCaSefAsxc9Jm4EJeQfHf2TYONJHcGFbuJxMMGsB92+0vIyt2ahy9whM2JghgMmOgWu1FxpC42zWn5uAiWG5uAKGpgk4Atbvxtj5c5G0XOF8EwctvQPD0RDZjmjAQ4Aoc7IrFC7ikRlBifuUNP+OIILQ+Bm7OR8snpj8sRQChUUUcVHpVjvn4LzN43deiBFEiT+\/fe6559wug4wC5ImwIOPhj4u8gOEcnkCOlpCA4XwdXBgjwg\/CIo1LACMIy2HKjok5+wt31GG+EAiLTrmFh\/bpJ3oHJgls2Mkxvkp+AIZa7GrBJobbW9UAAnEFDM1xf43Y9dlLH4\/luuLIABFQyj9cYiWzguhBpLBQI3aIB59d8eIBf0OQ+8JxJUcS0TeJysWgPzbx4\/R3YNhwsXniCAexj6BhEeXvpYXjI+r6bC\/HPMSY5zOyTmwoyCbRj8emHBbRoyx+TiaIFwIQK3AdfSEMiV+yMhyl+6Nv3x5z+M9\/\/tOyeaAugoI7R14sau\/AhObOW4\/+ZY3S+ZA9QXBW4+XoMwg7RCsbMYQvL05w1EamDXHjRViRAtYETAxrawWMfwsp+rYBAmTzzTdvOZ6JigSOfvhDapfzzlIBg2iAbCAdf6mMS7SQE0FK8ReDuehLCpnLqpyDIzpQ7OyeEFKkf0npcoeC4ndFoTsw1PXzh2i4XEY6EyKjLUSJv6znBQ3HI+wOGSvpT940gjiibx9wwRYSoi6lb9++TmBwnMbuhMJuCmzYeZFO\/u1vf+suN0NU7GIQJoyltHjy5syZezDcB2AnV+71VggLnDmmooAbhADBJ4FN6VhLFw8Ndq3BJobbW9U6ChgWHo5R8VdilGML4hN\/+\/e\/\/+2Oov0bfcQFWVkyhohd7oCxQYm+REC88CKAP+Lhvwhfjqwp\/k1F4pTFnZgiXohDf+xbLQbZJME9jA2hwAaFY\/DoW0iheyC85sxCy7EVx9XELVlpsoZkS8ACHCpd9q90hMT84Bc2SrSPEKINeAhRwrwQcqWvPSM6yLIyHjZ7xD0bPO4SIQThKi9g4FMECOP1n7RA5Pg3N0Nzj7oOfFx6hBTiZf8893442mfzSJ9kWxB92BC7wmnV3ujKaxCbgIlhWa2A8U1CFrwGTRaGP\/5uBUTA+asv3rFZuMlMEFzsUqLnzn4xJkBxenY71PWvT7PIk5nAofluCoSFo7NDgxxpj+MWyA4CZVFnN+Pf5NEs0j4t678DQ\/ucsfsdTqmAgey4YMc4uXeDgENwcS\/A79R8VqfUDIgNX4c58Qfi49wZcgIDdqCQE2fBEH9pgfR5jp0Ouysw4Z4Bwc\/dmXLfgYGwIArO6v13YJLAJiRgNNjRRq3YxHB7q1pHAUNT+DB+S0YEP0akwA8IFxYl7\/cccRAf3AFj44K\/lgoYvxjzLLtz+AHxH73HxsaFb8vgUxRin3t5\/ltMoRgkFhEKxAHiiDfgGLffYIQWcbIDiBiOtYhDNngIL4Qbc61VwCDEEATwgP80BJsO\/6kJNjl+jlETwo8IHjKXiB7inLuBZLX9kW70OzBkuNi0cbwLd0Uv+4bmHhIwIV72z5d+1wfRxvywCWNm4xu6vJzHQC6MgEEksMASjNxgj97g5ogAZ+DbBNw5IdjLfYAtKw7AJTWCnp2AFUPAEAgjwFEIAp57COy0ObIhM0hptvivNFs4jiwl340qfcsvjJDVMASaD4HCCBhScBwfsAtgBx8VMKTe2E2wM2eXT6aFb6tktbDrgqBIG1oxBAyBMAIcX3LsyS6VTAUbFjIf3AlrtvivJmC4mOvflAmjYjUMgeZGoDACxpuJIw1SmlEBw5kyF6O4Ke5v\/5OWJG1oxRAwBPKFAHfJuCTtMxUW\/\/myr82mOAgUXsD4D4ZxnurPdDmf5dXGSr\/jpzjuYTM1BPKHAK\/dksXkWMniP3\/2tRkVB4HCCxhe\/+O7BlEBw6u2XITzF7qK4w42U0Mg3whwsZsL8Vx+5vV6i\/9829tml28ECi9gMG8ohRz9zdH5dgebXV4R8L\/LJq\/z08yLOy+8+cLryLxt5ovFvwY9q9PMCOQ1\/k3A\/P\/P1fNGAq\/z+Uu8vK3gv0uC00J2\/pf9lXNiPqLEq89WxH1QyrD4X0\/IEhYhH8677\/K6Kd8u4cOHpb92oojxnyXf1Phes43X4l9j1dbXKYyA4euLfGuBj0dR+CQ9n9DnfX9eo2TR9a9R872WKMkVnfxb72bWQqMRKLoPE\/t8Nj76+2K4A8PFfYv\/Rnun9Z80AnmN\/8IImNY4SF6N3xpM7NnmQsB8uHZ7GXa1Y2dPZgOBvPqwCRiFf+XV+IqpW5WcIGA+XLshDbvasbMns4FAXn3YBIzCv\/JqfMXUrUpOEDAfrt2Qhl3t2NmT2UAgrz5sAkbhX81u\/Hq\/RZXXG+0KV2jaKs3uw40E3rBrJPrWdz0QyKsPm4BReEezG5\/x97tihmKmuip3D+5T9a0sXStWK00Emt2H08SqtC\/DrpHoW9\/1QCCvPmwCRuEdzW58EzAKI+e8SrP7cCPNY9g1En3rux4I5NWHTcAovKPZjW8CRmHknFdpdh9upHkMu0aib33XA4G8+rAJGIV3NLvxTcAojJzzKs3uw400j2HXSPSt73ogkFcfNgGj8I5mN74JGIWRc16l2X24keYx7BqJvvVdDwTy6sMmYBTeETJ+vd\/yUQwpdhW7xBsbslw9EPLhXE22zpMJYaeJf3tzr85GseZiIRDy4ViNZaiyCRiFMULGTyLDUW\/BUe\/2qv1uKAWkViVlBEI+nPJwmqq7EHaa+Lc395rK5LkbbMiHm3XCJmAUlgsZX0Ngim5aqkB29RYc9W7PBEwciza+bsiHGz\/C7I4ghJ0m\/k3AZNe+RRhZyIebFQMTMArLhYyvITBFNyZg4oBkdWMhEPLhWI0VrHIIO038m4ApmNNkbLohH87YcNXDMQGjgCpkfA2BKboxARMHJKsbC4GQD8dqrGCVQ9hp4t8ETMGcJmPTDflwxoarHo4JGAVUIeNrCEzRjQmYOCBZ3VgIhHw4VmMFqxzCThP\/JmAK5jQZm27IhzM2XPVwTMAooAoZX0Ngim5MwMQByerGQiDkw7EaK1jlEHaa+DcBUzCnydh0Qz6cseGqh2MCRgFVyPgaAlN0YwImDkhWNxYCIR+O1VjBKoew08S\/CZiCOU3Gphvy4YwNVz0cEzAKqELG1xCYohsTMHFAsrqxEAj5cKzGClY5hJ0m\/k3AFMxpMjbdkA9nbLjq4ZiAUUAVMr6GwBTdmICJA5LVjYVAyIdjNVawyiHsNPFvAqZgTpOx6YZ8OGPDVQ\/HBIwCqpDxNQSm6MYETByQrG4sBEI+HKuxglUOYaeJfxMwBXOajE035MMZG656OCZgFFCFjK8hMEU3JmDigGR1YyEQ8uFYjRWscgg7TfybgCmY02RsuiEfzthw1cMxAaOAKmR8DYEpujEBEwckqxsLgZAPx2qsYJVD2Gni3wRMwZwmY9MN+XDGhqsejgkYBVQh42sITNGNCZg4IFndWAiEfDhWYwWrHMJOE\/8mYArmNBmbbsiHMzZc9XBMwIjIkCFDZObMmbL22ms74E477TTp379\/C4gh42sITG0REbHfhRQHLaurQSDkw5o2mrnO7NmzZcSIETJnzhwZMGCADB48uGU67777rgwaNEjmz58v7du3lzFjxkjnzp3rGv8mYJrZe5p\/7HmNfxMwnwqYXXfdVQ488MCynhoyvgmY5g\/wvM8g5MN5n\/+CBQtk4cKFcu+998qyZctWEjDDhw+XNm3ayNChQ2X69OkyYcIEmTx5sgmYvDtFgeaX1\/g3AWMCJnYY224yNmQNfyCvBBYX2LFjx8rixYtXEjC9evWSSZMmSadOnWT58uXSs2dPl5Ht0KGDaz6EnWYDYzET11JWv54IhHy4nn2l2ZYJmE8FzEMPPSRrrrmmbLfddnLhhRdKx44d67oDi2NUO0KKg5bV1SCQVwLTzD1ap1TALF261B0XzZ0712VhKH379pWRI0dKjx49TMDEBdjqZxKBvMa\/CRgRIb282WabOce79tprZdasWXVPIcfxahMwcdCyuhoE8kpgmrlXEzBLliyRLl26rCRg+vTpI6NGjZLu3bubgIkLsNXPJAJ5jX8TMCXutmjRIpdCZke22mqrtRDYwIEDhT\/liiaFHMerTcDEQcvqhhAYN26c8Oell14KVc39z2s9Qmpt\/NsRUu5dK7MTzHP8m4ARkRdffFG23nprWbFihbvAN3XqVJk2bVqLQ4bUqwmYzMauDexTBEI+XBSgygmYYcOGSdu2bVsu8Y4fP16mTJlS1\/g3AVMUD8vmPPMa\/yZgROTYY491GZfVV1\/dpZMvvvhid3HPl5DxTcBkM2htVP+HQMiH847VvHnz3KcRPvroIzfVdu3ayXXXXSfdunUTXqMmw+Jfox49erR07dq1rvFvAibvHpbt+eU1\/k3AKPwuZHwTMAoQrUpDEQj5cEMHl\/HOQ9hp4t8ETMaNnPPhhXy4WadvAkZhuZDxNQSm6Kalit2BiYOW1dUgEPJhTRtFrRPCThP\/JmCK6j3ZmHfIh7MxyvijMAGjwCxkfA2BKboxARMHJKsbC4GQD8dqrGCVQ9hp4t8ETMGcJmPTDflwxoarHo4JGAVUIeNrCEzRjQmYOCBZ3VgIhHw4VmMFqxzCThP\/JmAK5jQZm27IhzM2XPVwTMAooAoZX0Ngim5MwMQByerGQiDkw7EaK1jlEHaa+DcBUzCnydh0Qz6cseGqh2MCRgFVyPgaAlN0YwImDkhWNxYCIR+O1VjBKoew08S\/CZiCOU3Gphvy4YwNVz0cEzAKqELG1xCYohsTMHFAsrqxEAj5cKzGClY5hJ0m\/k3AFMxpMjbdkA9nbLjq4ZiAUUAVMr6GwBTdmICJA5LVjYVAyIdjNVawyiHsNPFvAqZgTpOx6YZ8OGPDVQ\/HBIwCqpDxNQSm6MYETByQrG4sBEI+HKuxglUOYaeJfxMwBXOajE035MMZG656OCZgFFCFjK8hMEU3JmDigGR1YyEQ8uFYjRWscgg7TfybgCmY02RsuiEfzthw1cMxAaOAKmR8DYEpujEBEwckqxsLgZAPx2qsYJVD2GniHwETKr1795aJEyeGqtnPDYHYCIR8OHaDGXnABIzCECHjawhM0Y0JmDggWd1YCIR8OFZjBascwk4T\/5qva1uWpmCOleJ0Qz6c4lDq2pUJGAWcIeNrCEzRjQmYOCBZ3VgIhHw4VmMFqxzCThP\/JmAK5jQZm27IhzM2XPVwTMAooAoZX0Ngim5MwMQByerGQiDkw7EaK1jlEHaa+DcBUzCnydh0Qz6cseGqh2MCRgFVyDWvkqwAACAASURBVPgaAlN0YwImDkhWNxYCIR+O1VjBKoew08S\/CZiCOU3Gphvy4YwNVz0cEzAKqELG1xCYohsTMHFAsrqxEAj5cKzGClY5hJ0m\/k3AFMxpMjbdkA9nbLjq4ZiAUUAVMr6GwBTdmICJA5LVjYVAyIdjNVawyiHsNPFvAqZgTpOx6YZ8OGPDVQ\/HBIwCqpDxNQSm6MYETByQrG4sBEI+HKuxglUOYaeJfxMwBXOajE035MMZG656OCZgFFCFjK8hMEU3JmDigGR1YyEQ8uFYjRWscgg7TfybgCmY02RsuiEfzthw1cMxAaOAKmR8DYEpujEBEwckqxsLgZAPx2qsYJVD2Gni3wRMwZwmY9MN+XDGhqsejgkYBVQh42sITNGNCZg4IFndWAiEfDhWYwWrHMJOE\/8mYArmNBmbbsiHMzZc9XBMwCigChlfQ2CKbkzAxAHJ6sZCIOTDsRorWOUQdpr4NwFTMKfJ2HRDPpyx4aqHYwJGAVXI+BoCU3RjAiYOSIG6Rx99tDz++ON1a7Hev6cm7fGFfLhuQOWwoRB2mvg3AZNDx0hpShquCPFTyIdTmkrduzEBo4A0ZHwNgSm6MQETB6RA3SRs8tJLL9VthGmPL+TDdZtYDhsKYaexpQmYHDpGSlPS+lc1fgr5cEpTqXs3JmBE5N1335VBgwbJ\/PnzpX379jJmzBjp3LlzC9gh42scLI7lNGTX6PbquZjHmYu2bhI2qeec0x5fyIe1uOaxXhrxr4lp+2WOefSu1s9JwxUh38lr\/JuAEZHhw4dLmzZtZOjQoTJ9+nSZMGGCTJ482QRMhdgLBUvrQ7b1LWiCPk4v9Z5z2uPLK4HFsWGlumnEv1bA1GM+mjZCRw6aNqxOOghouCLET3mNfxMwItKrVy+ZNGmSdOrUSZYvXy49e\/aUmTNnSocOHZyHhoyvcbA4rq4hu0a3V89sRJy5aOsmYZN6zjnt8YV8WItrHuulEf+amE67Tj39OY9+kZU5abjCBExWrJXyOJYuXeqOi+bOneuyMJS+ffvKyJEjpUePHiZgytijWrCMGzdOBg4cmLIVV+1OE\/RxBhkiiHJtVcMi7fGZgClv7bTiP21x0u+KGVXduxZ\/jhMv9a6bFV6JM696jVnDFSF75jX+C5+BWbJkiXTp0mUlAdOnTx8ZNWqUdO\/e3fmr5hZ4HMe2uoZA2gjYkUF5xC3+0\/ZE668RCOQ1\/gsvYHCmUAq5EQ5nfRoChkA6CFj8p4Oz9WII1BsBEzAiMmzYMGnbtm3LJd7x48fLlClT6o21tWcIGAIZRMDiP4NGsSEZAgoETMB8+ho19zb8a9SjR4+Wrl27KuCzKoaAIdDsCPAatcV\/s1vRxl9EBEzAtMLqZGmuvPJKWbFihfTr10\/OPffcVrTW3I9+8MEHctppp8lzzz3n3tq66667mntCrRj9k08+KZdffrm89tprLrN38sknyzHHHNOKFu3RLCLQDPE\/ZMgQ90bl2muv7SAkRvv37+\/+\/\/nnn5dzzjlH3n\/\/fdl2223d96\/WXXfdhkB9wQUXyIMPPihvvvmmzJkzx8WNL5VwDn2\/J+mJVBrz66+\/LnvssYdsvPHGbgjM5eGHH3b\/3+gxJ41J2u2bgKkR8bffflv2339\/mTp1qnPUQw45RCCL3XffvcYWm\/uxxYsXy+zZs2XRokVy7bXXFlrAPPvss7LWWmu5LB5kdvDBB8uNN95oWb3mdvGVRt8s8Q8n7brrrnLggQeugj5vW5599tmy9957y0UXXSTrr7+++3sjyhNPPCFbbrml7LzzzisJmGo4h77fk\/Q8Ko2ZmD\/uuOPkgQceWGUIjR5z0pik3b4JmBoRnzhxojz11FMuA0Ph43fz5s1zO+8iF37\/0BVXXFFoAVNq\/yOPPFJOPPFE2WeffYrsGrmae7PEfyUB8\/LLL8uxxx4rjz76qLMLWY\/TTz\/dZUEaWcjeRjMw1XDOyuXr0jFXEzBZGXMjbVzPvk3A1Igm92Q4OvI7FtK0t956q9tpF7mYgFnZ+pDx8ccfL\/fdd59ssMEGRXaNXM29WeIfAfPQQw\/JmmuuKdttt51ceOGF0rFjR3nkkUfckZH\/4jjHSGQ\/XnjhhYbaqVQMVMKZFy1C3+9KayLlBMxee+0lm2++uayzzjouG8MRsuabQ2mNOS\/9mICp0ZIEFuWss85y\/50xY4bcfvvtcsMNN9TYYj4eMwHzf3Z844033C734osvLuzRYj68etVZNEv8L1iwQDbbbDM3AY52Z82a5UQLAoYPrd1xxx3uZ9zN2GWXXTIpYMrx7PXXXx\/8fldavlcqYDhOX7hwoWyxxRZCpos7R5dddpnwLZbQN8fSGnNe+jEBU6MlmyWFXOP0an7MBMz\/QsfZPbsu3m7hrpSVfCHQjPHP\/TR+TQpfHX\/llVecfz722GPOMHaEVLt\/lgqY0pb4qjv3i8444wz75ljtMJd90gRMjYD6y2V33323bLTRRnLooYfK4MGDC7\/TNgHzv7tZFocTTjjBXe62kj8EmiX+X3zxRdl6663dcTf39HjpYNq0ac4gfHGct5D8JV7eQGr0m5SlYqAazln5fk\/pmMl6cVwMnrxVBRecd9557g5cVsacl4g0AdMKS955550yduxY9wsgDzroIPcWUpHLvvvu61Kn7PQ23HBD9\/rwgAEDCgfJNddc49Lz\/jVKAGBhMDGTL1dohvjnCJOMy+qrr+6OLzjOZMGl8LYcmy4+gbDNNtu4FxLIFDSi+Ls6b731ljvy2mmnneTqq692Q6mEc6O\/31NpzNyHvPTSS4VfU8Er1PwqmlNOOcXNpdFjboRtk+zTBEyS6FrbhoAhYAgYAoaAIZAIAiZgEoHVGjUEDAFDwBAwBAyBJBEwAZMkuta2IWAIGAKGgCFgCCSCgAmYRGC1Rg0BQ8AQMAQMAUMgSQRMwCSJrrVtCBgChoAhYAgYAokgYAImEVitUUPAEDAEDAFDwBBIEgETMEmia20bAoaAIWAIGAKGQCIImIBJBFZr1BAwBAwBQ8AQMASSRMAETJLoWtuGgCFgCBgChoAhkAgCJmASgdUaNQQMAUPAEDAEDIEkETABkyS61rYhYAgYAoaAIWAIJIKACZhEYLVGDQFDwBAwBAwBQyBJBEzAJImutW0IGAKGgCFgCBgCiSBgAiYRWK1RQ8AQMAQMAUPAEEgSARMwSaJrbRsChoAhYAgYAoZAIgiYgEkEVmvUEDAEDAFDwBAwBJJEwARMkuha24aAIWAIGAKGgCGQCAImYBKAtVu3bvLhhx+u1PLaa68tf\/\/7392\/XX311TJmzBj3\/+ecc46ceuqpLXVXrFghP\/\/5z+W2226Tt956S9Zbbz3Zaqut5Ac\/+IHstttuUq5t\/\/DnP\/959+wBBxzQ0t5aa60lX\/jCF+Skk06Sww47bKUx3XzzzfLrX\/9aHnjggSAKhx9+uDz99NMr1dt\/\/\/3dXLTl7LPPlscff1xef\/1198if\/\/xn2XTTTVsenz59ulx77bXy0ksvySabbCLf\/\/73ZcCAAe7ns2fPlsGDB8u\/\/vUv2XHHHR1+\/tnjjz9eunTpIkOHDq04lBCuPBjFtl27drLZZpvJ7rvvLqeddpr7\/ziF8V5xxRXyl7\/8RdZcc03Zeuut5brrrnPtYKNp06bJa6+9Jquttppss802MmTIENc\/47zgggvcz\/GZH\/3oR3LUUUe1YMBc77\/\/\/tjjiTP2LNYldqJ+7ce41157yfjx41sVVxtttFHZtn0fRxxxhPTv318dV8TqBhtsICNHjgxCGYqJOXPmyPDhw53\/r7\/++vLtb39bzj33XOdTmrJ06VK55ppr5M4775R33nnHccn5558vX\/va16o+Tv1f\/vKX8s9\/\/lM+85nPOO6hnXqWanP7xz\/+4eLnr3\/9q7z33nuOD5j7oEGDZI011nBz+eEPfyhPPfWUfO5zn5Mf\/\/jH0qNHDzc84uvee++Vu+66y8VXuZImH3zlK1+R\/\/73v6sM45JLLpFjjjlmlX9\/\/\/335ZRTTpEXX3zRrSMbb7yx7L333s5u2KKeZcaMGTJ69Gh59dVX5bOf\/aycfvrpLesEPPOrX\/3KjeOjjz6SLbfc0q1VBx10kBsCa9SVV14pn3zyiRx44IFy6aWXOrw\/+OAD2WeffeSyyy5z\/02ymIBJAF2\/EHbu3NkJEAqOd9NNN7n\/33fffYUApeDc99xzT8sobr\/9djnvvPPc4tWzZ09ZtGiRzJ071y1kJ554orCA\/ec\/\/3EB8dxzz7nnunfv7oK6Q4cOzsE80UPub7zxhrzwwguu3syZM4Ux+UJ7iJsLL7wwiIIXMCzovvTu3Xsl8RVq5Bvf+IZ86UtfkgcffNAt1FEB89hjj7lghqhYsO+77z6ZN2+eI82+ffvKoYceKuuss44Ltm9+85uuDoKGIAMvRBgEX6mEcOW5qN2WLVsmr7zyivBfhNKUKVOkY8eOoSm6n7PYQrYsHn369HEExL9BspAA84GY6A+CRshsuOGG8vDDD7uFChsjLp988kk3f\/6NudMmohExWrQSFTA77LBDy\/T5f8Rfa+IKkkUUUP797387u4N3165d3b+xeOyxxx6quMJfGNNPf\/pTFXlXi4nFixfL17\/+dXn77bedkEbYL1iwwC0yZ511lsoFiBc2Bdtuu63sueeejoNol3ip5M+33HKL44R1113X+erqq6\/u+IoYilvw1Z122slxV7SE5vb73\/9ezjjjDIclmwm\/yfIbPvAlRuA0bEc8IVjefPNN2W+\/\/eTGG29sETTlxpwmHzD3qIB55JFHHP8hEMst8Gzwjj32WLexWb58ueNLnkfUeD\/V2oH1o1+\/fnL55Zc7O0SLjyk2ufgGHEx9cEXgjhgxwvFer169nA8+++yzbtyIW8b21a9+1fkiP2e8bIaJE4QZvsrfky4mYBJA2C+EEydOFBb5aGHXwSLEooZaZSfBboEdOoUdGU4TVeeoX3YhUcJh0YPcKCyCkI1fPL2AwYlwOOpBfGRL6JtCQCCQ2K1AjmQ2UNDPPPOMIywcERLzogABg2DyYqgUNrIlf\/jDH2TcuHFOjVcrfkcSFTDDhg2T3\/zmN3LmmWe6XRbCBPLbfvvt5e6773ZByPyvuuoqNx+EAcEMWQ0cOFAOOeSQqn1qcC2128svv+x23mCHeGAxIIj5f4RIpczVySef7AQYBEBglxbmTdB7+++8887OTpMnT5Z3331XTjjhBGFnxO6SbAz2JTt1\/fXXO8Ju06ZNAl6b7SajAga\/Li31iivIGb9iwScL5ktp\/5Xiiizl0Ucf7bKVxCSLAostu1hEKsKbhYgYi5ZyMeEXWTgELkHQf+tb33Lt4hssPKG4Iz7mz5\/vnqcddsXsqvGxchlLRDf1Fi5c6HyQzGZpYZFj543IYLdN3Fx00UUtHBatjwCCS9hsREtobv\/zP\/\/jMPKcRxYKUQJX8V\/mgWABXzZ3+AQxQlaG7NRPfvKTzPBBdCD4BVxKtoMNCxvPUGEtQAxEM95av8I+bHC9uIj2BW5Tp051mCEWEbcXX3yx49cbbrjB8T0bTgQk\/k7\/8D8i8rjjjnPcTH2yc2zKEWTwIvyPLfj\/pIsJmAQQ9gshxvcZGNQqQUhalKME0tKUSZMmrbSjGjt2rFukcXAcgR3ILrvs4pwoWrQChsUXx8ORoxkYdgHsDDjigLTIaqD8GRfiBseGeLyKJuggTQQNCy9qHgL0Thoi0hBZe2JlR4II+tnPfuZ2jpA0QUOwE3Aff\/yxIzUCDAJDRCD4KqWKfb8aXMsJT9KkzJN5g5VGwGBrSB6y5Rl8ACHDwlVayJCx24Gs\/\/SnPwnHGRxZ\/fGPf3RVOY4gI0MWgIUQmxSxRAUEdvKFrB1Ho\/WKK62AqRRXiFwyZuxiEVWk28mmQvhk1FjAyBiVZtGqifqo2PD18HtiLxR3CH8WIrKU3\/nOd1zG9NFHH20RApWEID5LXzz7xS9+0Y3Z+x7+CZfsuuuuTugRi9RnA4O4YrdOhoXy3e9+VxDo\/picOOKP37BUm1t0bPQPVxJHbAzAn+NleLBt27ZCfNP39773Pbf5iR5Nl4uXNPkg2j9CgcW99OpA6RjJvHBMznESGzgyewhPxKXGr8CFQmYK\/if76zfTYANmbALJ7tMu\/ILfkunm52yyooXxIFLY1PmsOBzNWkVB9MDX+CM+USpYk+IsEzAJIFvungqLFKSGWkUo4DQ4BUQGQRD8FIQGzs0O3hcCE4LG4XzRCBhfFxGAeDryyCNbnkcwkBbm\/gABRWCRFcApKRAv43zooYfcMROZCEQCRzzcY2HhhUB9RgBBRPoWsuIuTrVSjqzJANEn86eg+L0oQQRAevwbAcn\/0z+7SzAlg8UfAg9i9Kn\/6Bg0uJYTMOx0fBaFXTTtsHhgE3bDpQWi8dk0jvSog9CCiEqzU2TWIFwWNsQZYsUXhBo7SbItEDZn1NTBjiwQ7Gy1xwgJuHjqTVa6A+Nxq1dcaQRMtbhi04Efszkge0bmgwUfOxMfZEe32GILJyKipVxM+MXOZyWpzz0PfIkMBpuIUNyRRcGviJ1oXBEjxH1poT1inUI8kIFBOMAhCAMWPmKcIzY4ir+zC\/\/d734no0aNchxT7r6c78fPRTM3\/wyZHsQ\/98fIikXFCfGIeGJ+4E521N8z43mOOHzWuRF8EO0T+7PQ+80KGblKhU1l9LgfbibrxcY25Felz5b24bNxtMkJgP87YoY1hqwQXBct3OdirSDOEKw+g0hfS5YscVcksBNZezaUHFmROcaHyOxXO95vDZmYgGkNehWerXSE5NOHBD\/\/j4DhGIfdCotcdGdJChehAHmwGyco+bsvGgEDmbIAc67Pghe9hMfiz8IMWf3iF79wCr1c8alHxuqdlgUUNQ9pQILlBEM1WMuRNfUJJu4DQUqcsbI7YwFnx1GacocAyUpxdsvu7tZbb3UBBklEU\/+l46iGqyYDo3EXUrbMwRM6uEM+Bx98sPuvnyvzI8g5cqq0Y+HYgOMxcIHEEYdkZLBdpTN0zRibrU61I6R6xpVGwFSKK+KCRQExzeJD9sFn0krx\/tvf\/rbShcxqGRg2Of6eD4KCBcNnYDR2ZGEig8mGgJhlcWEh8nfyom34ozj+jbsXnTp1cr4Gj+DP9F\/puNb7Mbt3RBYF8b3ddtu1PMMlYv74DExobgg1MqCIPo6YKx1LcG8HrmJO3P3jWATOYtPgM5vlsEqDD3y\/\/j4SWW5soCmMj80mosAf7YT8ioyf3wBzXxKRD1+SGaHwEgSZXp+BgTfBrFwGhg0Zl4dZh9iEc4xd7iIxaxjZczZVZKnpn7pwHAKy2gsWGhwq1TEB0xr0YgoYH7TlHiP1hqNAFGRk\/NmvV8XRt5h4XiNgOBcmi8JFQRzMixHSrxxv+OwKCz73SDiu4r++cJTEcQhCgUDyY0Jo+Hsb\/ljKv13Ev4cuu5Yja8gGgvVZFxbnSkQLQUJ8ZK0gNeZFloaA4ZgFzEqLBlfNHRjuqNAvO4pKN+y55\/DEE0+4rBm7UbJajIsUPjsZbIdwJB3LbgrbVyosHuy6IaEvf\/nL7q4BAobzZxaMaNYmAVfOTJPVBEw940ojYCrFFZsQ7OyzK8QzwhobsmGIxhWLUfRNonIx4Y8vffbW34Eh44B\/kf0IxR07ZN8P4oUFBf\/lThiZIcQQ8U8dMhhkBVng4AtED4IZXyUzgxhg40KMk4HkyCB6dEtmpPTOTKU7MJq5+bc18Xt2\/WQfyhXEEkcgjI9MEfHHxorFlzndcccdjtuiJU0+oF\/w5MgNHi3d9HFMR8YDoYfwZfMTzVj4+0KISURlHL+qdgcGkcfxlM9i+jsw\/p4RWWAyZohl7MjmC3zLFWw1a9YsJ3TgOd7CY2NMJpK5lxPL9SAXEzD1QLGkjXI7eYKJ4CcbQuoXA1MIPo4QyLAQVKTgID12GgQsOyKcPvq6aBwBQ11\/AY7dIcRByhDV7S+hcl+Dow7EDilg+uX1SX7O2Hidm\/Nvxs1OCAFBXUgX8iPlGDqLZxxcRmO+7CYoBAXCjJ0dl5QhHu77+AuEZF8IXv96JM8gdEjDM04W8t\/+9rfu0i+kgFBAvJDKLC0+6KvhGr27BNnzNgr\/hZhZnMBFcwfGX0Dm9U7uIDAHduf8FxuQiSHzwkKEEPEFQokSLel\/0q+0x64HIYqghBTY1ZDZKZceT8ClG95kJQFT77jSCphycQXZs7CzY6Y8\/\/zzzv7YGd\/Gn\/FPdsUcfYZiAvGB0IEzWKBZ5DhKjL6FFIo74pO+iFXGwy6bGMa3WCT9AufvmjEmfx8NX2XXzuLTvn17F1cce\/hL6mRx4ATGR6aJbAv\/Fi2VBAyLY7W5sVAjtiiMwfMlm7vS19M5wiLzyoIJL7EAMxZ8A85BUJZ+BiFNPmAO3q88B0cxYqwINH8vhhcrqM+8\/VtIvKrMPSq4XONXvv1qAoZ2yCbCwXALAoQjeuyNHb0f0Bai1V84Zq2IvpzAhph\/g9\/InJNtQVjjd\/g9axfjTqKYgEkA1XICxp8tE4gYl90ThUBG2OA4CAscFsOzwJHpgPxwJpwievarzcDQBw6GE0GIvOlCgJAlYefqC9kAMh68hcRY+Dn9Qg4QLs6MwIIgmAM\/Yxfnd0UhIqUf5okYKi2INLI97PQgaXaFLNSkI0vf4gIb\/iBYCChIivQ6mSBIqtzrgvTH2EO4Ru8ucczHHRZ2TWQ5PPYaAUN\/pLQhVOZLypwUrr\/D5N8MKcWB+n6nzq4FUuFoyb\/VhXDED2gTEQeRlx6tJeDOmWiykoCpd1zFETDRuGKXz4KLGI\/ejWJ8XMb0dwoQ0Bwh+O9\/VIsJOILYIAbxO8QRggh\/92+iheKODQj1GSvteS7xWdJyAgbfI97ZHMAZxAULvj8q5hkuwSJo4CiOI8ja4Kul2ddKAganqjY3fy+v1PlKPzuBIOTiKfHPhoHCAuy\/T8Xmpty3VtLmA+IVHyZrFb3LyHhLBQwZMo6b\/Gcc4FhwhIf8uhHyK42AoQ5H0\/47MPAd4tjflayU2eQIi0ywL6XfPUJMgjs+y3E6F31Dl6prJZnCCBh2HjgKQcPt9eidA78AokQppFb9RbZagc3qc4gD7t2wey\/qGy1ZtY2NKzkEqsU\/x4IQLq8bk2lAcEQvUGpGRfsQPxsA\/0kDzXNWxxAwBGpHoDAChh0IRxNcsGPXXipg2GmHvl9SO8zZeZKjDLIDpIHr\/VXH7MzSRmIIrIxAtfgnvU1Gg+wWO\/8JEya4TGWcgoDhYm657\/7EacfqGgKGgB6BwggYDwmpT1KkRRUwetewmoZA\/hAoF\/\/cS+DyIZck\/ZuBHEmQUrdiCBgC2UXABIyIOyPmjRxu4nMTnIu0oTdpsmtSG5khYAhUQqBUwPhvZnCXwt8r4b4B94uil8cNUUPAEMgeAiZgPr3k6m+pc8mK29hxU8jZM62NyBAwBEoRKBUw3Anj1d+ogOGSJZdxuYBoxRAwBLKLgAmYEtvwSjGXXCE0\/40DfrdJ9CNy2TWnjcwQKI+A\/306RcenliMki\/+ie03zzz+v8W8CRsS94sjn3\/nmBxf4+D1A0a+58juNyv0COe\/W\/LzfFTOCXn734D5V2wk2kHAFPgQV\/ZBdwt0l2nxe5lKveYR8OFFjZKjxcgKG10V5PdVf4uUbSf5bRQw9hF29bBSXT7SwGu9okapcr942bv2IVm4hNL6QD9d7PGm1VxgBw1cseTWaL01S+Losv1SRbxzw5gAZF76pQTqZ7y5g8CihFEHApOV01k\/6COSVwLRIVot\/XqNGuPvXqPkuRvTXY6SNnXZDpJ171gWMdh5Wr3YE0vbh2kca78nCCJh4sKxcO2R8LeEYkbTGCvZsaxAI+XBr2s77s2ljp+UTLe7GO1qk8lsvbR9OC0kTMAqkQ8bXEo4RiQJsq5IIAiEfTqTTnDSaNnZaPtHCa7yjRSq\/9dL24bSQNAGjQDpkfC3hGJEowLYqiSAQ8uFEOs1Jo2ljp+UTLbzGO1qk8lsvbR9OC0kTMAqkQ8bXEo4RiQJsq5IIAiEfTqTTnDSaNnZaPtHCa7yjRSq\/9dL24bSQNAGjQDpkfC3hGJEowLYqiSAQ8uFEOs1Jo2ljp+UTLbzGO1qk8lsvbR9OC0kTMAqkQ8bXEo4RiQJsq5IIAiEfTqTTnDSaNnZaPtHCa7yjRSq\/9dL24bSQNAGjQDpkfC3hGJEowLYqiSAQ8uFEOs1Jo2ljp+UTLbzGO1qk8lsvbR9OC0kTMAqkQ8bXEo4RiQJsq5IIAiEfTqTTnDSaNnZaPtHCa7yjRSq\/9dL24bSQNAGjQDpkfC3hGJEowLYqiSAQ8uFEOs1Jo2ljp+UTLbzGO1qk8lsvbR9OC0kTMAqkQ8bXEo4RiQJsq5IIAiEfTqTTnDSaNnZaPtHCa7yjRSq\/9dL24bSQNAGjQDpkfC3hGJEowLYqiSAQ8uFEOs1Jo2ljp+UTLbzGO1qk8lsvbR9OC0kTMAqkQ8bXEo4RiQJsq5IIAiEfTqTTnDSaNnZaPtHCa7yjRSq\/9dL24bSQNAGjQDpkfC3hGJEowLYqiSAQ8uFEOs1Jo2ljp+UTLbzGO1qk8lsvbR9OC8nCCJjZs2fLiBEjZM6cOTJgwAAZPHhwC8b8NtpBgwa1\/DbaMWPGSOfOnVt+HjK+lnCMSNJya+unFIGQDxtilRFIGzstn2htZryjRSq\/9dL24bSQLIyAWbBggSxcuFDuvfdeWbZs2UoCZvjw4dKmTRsZOnSoTJ8+XSZMmCCTJ082AZOWF1o\/iSOQVwJLHDgRSRs7EzBpWLVYfaTtw2mhWxgB4wEdO3asLF68eCUB06tXL5k0aZJ06tRJli9fLj179pSZM2dKhw4d3GMh42sJx3ZCabm19WMZmPr5QCj+69fT\/7ak5RNtv8Y7WqTyWy9tH04LycILmKVLl7rjorlz57osDKVv374ycuRI6dGjhwmYtDzRy60JVQAAEglJREFU+kkUgbwSWKKgfdp42tiZgEnDqsXqI20fTgvdwguYJUuWSJcuXVYSMH369JFRo0ZJ9+7dTcCk5YnWT6II5JXAEgXNBEwa8FofKSCQ1\/gvvIDBdzRHSAMHDhT+lCvaHZOlclOIVOtiFQTGjRsn\/HnppZcMnRoQSJv8tXyinYrxjhap\/NZL24fTQtIEjIgMGzZM2rZt23KJd\/z48TJlypQWG4SMryUcI5K03Nr6KUUg5MOGWGUE0sZOyydamxnvaJHKb720fTgtJAsjYObNmyf9+\/eXjz76yGHbrl07ue6666Rbt27Ca9RkV+bPny\/t27eX0aNHS9euXU3ApOWF1k\/iCOSVwBIHTnGJv95jMAFTb0StvbzGf2EETGtcOGR8LeHYTqg1VrBnW4NAyIdb03ben00bOy2faHE33tEild96aftwWkiagFEgHTK+lnCMSBRgW5VEEAj5cCKd5qTRtLHT8okWXninnqV3794yceLEejZpbSWMQNo+nPB0Wpo3AaNAOmR8LeGYgFGAbVUSQSDkw4l02iSNDhkyxH33ae2113YjPu2009xxsy9pY6flEy288E6\/K2ZoqwfrGY8FIcpchbR9OC0ATMAokA4ZX0s4FvgKsK1KIgiEfDiRTpukUQTMrrvuKgceeGDZEaeNnZZPtPCagNEild96aftwWkiagFEgHTK+lnBMwCjAtiqJIBDy4UQ6bZJGTcDEM5TxWDy8slA7r\/FvAkbhXSHjxxEwiu7Ezpg1KFmdOAiEfDhOW3mri4B56KGHZM0115TttttOLrzwQunYsWPLNNPGTssnWjtYBkaLVH7rpe3DaSFpAkaBdMj4WsLREontcBRGsSqxEAj5cKzGclaZX\/S62WabuVlde+21MmvWrFi\/zLXecGj5RNuvlnfitGcfRdSilY16eY1\/EzAK\/woZX0s4WiIxAaMwilWJhUDIh2M1luPKixYtcr\/Mld+Nttpqq7mZgl21L3HXGw4tn2j71fJOnPZMwGjRany9PH+J2wSMwr9C5K8lHC2RmIBRGMWqxEIg5MOxGstZ5RdffFG23nprWbFihUyYMEGmTp0q06ZNa5ll2thp+URrBi3vxGnPBIwWrWzUS9uH05q1CRgF0iHjawlHSyQmYBRGsSqxEAj5cKzGclb52GOPdRmX1Vdf3f1i14svvthlXXxJGzstn2jNoOWdOO2ZgNGilY16aftwWrM2AaNAOmR8LeFoicQEjMIoViUWAiEfjtVYwSqnjZ2WT7Rm0PJOnPZMwGjRyka9tH04rVmbgFEgHTK+lnC0RGICRmEUqxILgZAPx2qsYJXTxk7LJ1ozaHknTnsmYLRoZaNe2j6c1qxNwCiQDhlfSzhaIjEBozCKVYmFQMiHYzVWsMppY6flE60ZtLwTpz0TMFq0slEvbR9Oa9YmYESktZ8S1xKOlkhMwKTl\/sXpJ68EloYFQ9gdffTR8vjjj9d1KPX+9H+926vnZO27V\/VEs3xbIR9OfgTJ9GAC5lMB05pPiZuAScY5rdX6IZBXAqsfQpVbCmGnjX\/tWLUbnTy1ZxkdrTVrqxfy4dpabfxTJmBMwDTeC20EiSOQVwJLHLhPvwNTbYE1AdM6K1jGuXX4aZ7Oa\/ybgPlUwLTmU+JaAtPurCygNSFpdeIgkFcCi4NBrXVD2GnjX9u\/lify1J5lYLTWrK1eyIdra7XxT5mAEZHWfkpcS2BaYjIB0\/jAyNsI8kpgadgphJ02\/rVj1fJEntozAaO1Zm31Qj5cW6uNf8oETIkNavmUuJbAtMRkAqbxgZGnEeT5U+Jp2ClE\/tr4145VyxN5as8EjNaatdUL+XBtrTb+KRMwItLaT4lrCUxLTCZgGh8YeRtBXgksDTuFsNPGv3asWp7IU3smYLTWrK1eyIdra7XxT5mAEZHWfkpcS2BaYjIB0\/jAyNsI8kpgadgphJ02\/rVj1fJEntozAaO1Zm31Qj5cW6uNf8oEjMIGIeNrCUxLTCZgFEaxKrEQCPlwrMYKVjmEnTb+tbBpeSJP7ZmA0VqztnohH66t1cY\/ZQJGYYOQ8bUEpiUmEzAKo1iVWAiEfDhWYwWrHMJOG\/9a2LQ8kaf2TMBorVlbvZAP19Zq458yAaOwQcj4WgLTEpMJGIVRrEosBEI+HKuxglUOYaeNfy1sWp7IU3smYLTWrK1eyIdra7XxT5mAUdggZHwtgWmJyQSMwihWJRYCIR+O1VjBKoew08a\/FjYtT+SpPRMwWmvWVi\/kw7W12vinTMAobBAyvpbAtMRkAkZhFKsSC4GQD8dqrGCVQ9hp418Lm5Yn8tSeCRitNWurF\/Lh2lpt\/FMmYBQ2CBlfS2BaYjIBozCKVYmFQMiHYzVWsMoh7LTxr4VNyxN5as8EjNaatdUL+XBtrTb+KRMwChuEjK8lMC0xmYBRGMWqxEIg5MOxGitY5RB22vjXwqbliTy1ZwJGa83a6oV8uLZWG\/+UCRiFDULG1xKYlphMwCiMYlViIRDy4ViNFaxyCDtt\/Gth0\/JEntozAaO1Zm31Qj5cW6uNf8oEjMIGIeNrCUxLTCZgFEaxKrEQCPlwrMYKVjmEnTb+tbBpeSJP7ZmA0VqztnohH66t1cY\/ZQJGYYOQ8bUEpiUmEzAKo1iVWAiEfDhWYwWrHMJOG\/9a2LQ8kaf2TMBorVlbvZAP19Zq458yAaOwQcj4WgLTEpMJGIVRrEosBEI+HKuxglUOYaeNfy1sWp7IU3smYLTWrK1eyIdra7XxT5mAUdggZHwtgWmJyQSMwihWJRYCIR+O1VjBKoew08a\/FjYtT+SpPe1cGlWvd+\/eMnHixLp1f\/TRR8vjjz9et\/ZC4wv5cN0GknJDJmAUgIeMryUwLTGZgFEYxarEQiDkw7EaK1jlEHba+NfCpuUJa68yAklgWM8sURI+U218IR\/W+lLW6pmAEZF3331XBg0aJPPnz5f27dvLmDFjpHPnzi22Chlf64zaoKJePUtIndezr2hb2l1Go8aX1Lzr2W69MAz5cD3H3GxtpRX\/Wly0PGHtmYDxCIQ2vXmNfxMwIjJ8+HBp06aNDB06VKZPny4TJkyQyZMnN1TA9LtiRpCftEQXcu5gRzVWiCPs6rm7qXG4mXysXhjmlcDqYbS04l87Vm1cW3smYEzArFixQhsIea3Xq1cvmTRpknTq1EmWL18uPXv2lJkzZ0qHDh3clEPkH2eRqbcw0bbXCIEQB5dGjK8Z\/LleGIZ8uBmwSGqMacW\/dvwmYLRImYAxAVNwAbN06VJ3XDR37lyXhaH07dtXRo4cKT169CiUgBk3bpwMHDiw9ezxaQv1WnxrGVC951LLGOrxTL0wNAFT3hppxr\/WH0zAaJEyAWMCpuACZsmSJdKlS5eVBEyfPn1k1KhR0r17d+cf2nsIrQ87a8EQSAYBu2dUHleL\/2T8zVrNFgJ5jX+7AyMioRRytlzRRmMIGAL1RMDiv55oWluGQHoImIARkWHDhknbtm1bLvGOHz9epkyZkp4VrCdDwBBoGAIW\/w2D3jo2BFqFgAmYT1+j5u6Hf4169OjR0rVr11YBaw8bAoZAcyDAa9QW\/81hKxulIRBFwARMK\/yBLM2VV14pXCPq16+fnHvuua1oLflHL7jgAnnwwQflzTfflDlz5risky+V5hL6Rkbyo165h7fffttlyp577jlZY401ZM8995SLL7645QJ2s8zDz+qoo46SV199VZYtWybbb7+9u3u1ySabuB8321zS9oVG99eo+J89e7aMGDHCxfCAAQNk8ODBLVBUi9fnn39ezjnnHHn\/\/fdl2223dd+7Wnfddav6Wi0Y33bbbUIW+4MPPpCNN97YZbh33nln11QWxsenMtikfvLJJ7L22mvLaaedJocddlhmxhfF\/OWXXxbuZJ544oly1llnZXKMtfhIvZ4xAVMjkiyk+++\/v0ydOtUF6SGHHCJDhgyR3XffvcYWk3\/siSeekC233NKRSVTAVJtL6BsZyY96VQEDEe+xxx6yePFiF9j77ruvHH\/88dJM8\/CzWrBggXTs2NGJ4EsuucQJGTBvxrmk7QuN7K+R8Y\/PLFy4UO69917nL1EBUy1eebvy7LPPlr333lsuuugiWX\/99d3f6z0XNkmI8U033VRmzZrlPhL65JNPyuqrr171m1tpje+dd96R9dZbz23gwJJ+H3jgAbdxyAJ+3q\/hhO9+97turPC2FzBZGmMjY5C+TcDUaAF+L8ZTTz3lMjAUPn43b948ufzyy2tsMb3HeKU2KmCqzSXrFxzZSSFkyMo08zxYiBAwkBYE1cxzSc+TG9dTFuJ\/7NixzvejAqZSvH788cdy7LHHyqOPPupAI\/5PP\/10l5FNci7\/\/e9\/ZbvttpNnnnnGZXuyNr5XXnnFZc\/vv\/9+J7iyNL6bb77ZZazAkOIFTJbG2LgI\/N+eTcDUaAEWThYbdjAUPnx36623yo033lhji+k9VipgKs2FNHDoGznpjXrVnj766CM54IAD3LELrwk26zwOP\/xwIasE1jfddJNssMEGTTuXRvpDmn1nIf5LBUy1b9oQKxwZ+S+Mc4xEJvaFF16o6Gv14LJbbrlFpk2bJhwrZWl8HCOxUXjvvffcN7+IwSyNj8zQmWee6XC76qqrWgRMlsaYZrxV6ssETI1WgMCiqnjGjBly++23yw033FBji+k9Vk7AlJvL9ddfH\/xGTnqjXrknvt9xyimnuFQ1KWpKJZtkeR5+VuykudfAcRJn8s08l0b5RJr9ZiH+SwVMtW\/aLFq0SPi44x133OFgYme\/yy67tAiYJLjssccek\/PPP1\/IJODXWRsfcyYTxQVuuJsMUaVvgqWN3wknnCCnnnqq7LDDDitxQRYxTDPuSvsyAVMj+kmmXWsckvqxZj9C4riF3cnmm2\/uzvJ9afZjFy5nQvj33HOPHSGpvbkxFbMQ\/3GPkI455hhBVPiFO8kjJO68cOTBhg6+8aXa8Uea44t6Dffn6Hu\/\/fareoSU5vh23HFHd8GY8uGHH7r\/Hnnkke6oPIsYNiYK7QipZtz9xbe7775bNtpoIzn00EPdWXSWL\/H6yZYKmGpzydo3MvhdVRAjwX3ZZZfJaqut1mLDZpoHg+YiJql8Luixs7r00kvdmxuk+pttLjUHUpM+mIX4LydgqsUrb7PwFpK\/xEvGgTcn6z2XZ5991t2vIfNZ+jmKLIyPNxjJtKy11lru3iLCgKM1eDEL4ysNidJsXxbH2KgwtgxMK5C\/8847BRJhUT3ooIPcW0hZLozvoYcekrfeeks222wz2WmnneTqq692Q640l6x9I+Ppp59259W8McBbDRR2Thy\/NNM8GCuvT3MMxgKy5pprunQx5ORfo24Wm2TZ55McW6Pin0W3f\/\/+wr0WSrt27eS6666Tbt26uaOhSt+0QViwyUIkb7PNNu4FBN5EqhY3teDH2HjBwbdNG9zl4JflZmF8cDbjoTBGjmwPPvhg9\/csjC8kYLI4xlr8pB7PmICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4ImICpB4rWhiFgCBgChoAhYAikioAJmFThts4MAUPAEDAEDAFDoB4I\/D814Rw2uj\/8+AAAAABJRU5ErkJggg==","height":337,"width":560}}
%---
