%[text] %[text:anchor:T_360F0929] # Artifacts in resting state EEG data
%[text] This script uses ArtSim's option to load data (from a Fieldtrip or EEGLAB data structure) and investigates whether (simulated) tACS can be removed and the true underlying signals recovered. 
%[text] Specifically, I used 4 minutes of (eyes open) resting state activity recorded at electrodes near Oz.  In these data the 1/f decay rate as well as the strength of the alpha rhythm varies considerably across participants. 
%[text] The simulated tACS signal applies a 1 mA current at 10 Hz  that results in a 15 mV artifact in the EEG signal.
%[text] I then apply the three artifiact removal algorithms and compute the residual artifact (expressed as the percent error in the estimated power) near 10 Hz (directly affected by tACS) and near 5 Hz (control frequency).
%[text] %[text:anchor:H_D32E8D0B] ## Setup Simulation
%[text] First define the stimulation and recording parameters.
run('../startup'); % Make sure all toolboxes are on the path
defaults = [];
defaults.reproducible              = true;    % Fix the random number generator to get reproducible results
%[text] Stimulation and Recording Equipment:
defaults.stimulatorSamplingRate    = 1e3;     % 1 kHz 
defaults.currentResolution         = 1e-6;    % 1 muA
defaults.highPass                  = 0.01;    % Analog filter of the recording amplifier
defaults.lowPass                   = 100;
defaults.recordingSamplingRate     = 250;
%[text] Experiment:
defaults.duration                   = 4*60;   % 4 min of data
defaults.tacsAmplitude              = 1e-3;    % 1 mA 
defaults.resistance                 = 15;     % 15 mV artifact per mA current
defaults.tacsFrequency              = 10;      % 10 Hz tACS
%[text] And setup the base parameters for artifiact removal
ar = [];
ar.recordingSamplingRate = defaults.recordingSamplingRate; % Match to the recording
ar.segmentDuration       = 3;
ar.tacsFrequency         = defaults.tacsFrequency;

%[text] Setup the parameters for three artifact removal algorithms to explore.

fastr = [];
fastr.name = "FASTR";
fastr.options.slack = 5e-3;
fastr.options.nrSegsPerWindow = 10;
fastr.options.pcaNr = 2;
fastr.options.pcaVarExplained = Inf;
fastr.options.pcaNrSegsPerWindow = Inf;
fastr.options.referenceSegment = 1;

fbar= [];
fbar.name = "FBAR";
fbar.options.frequencyLimits = [-1 1];

rbar = [];
rbar.name = "RBAR";
rbar.options = [];
%%
%[text] %[text:anchor:H_DDD99D7B] ## Estimate Alpha Power
%[text] Create the artSim object using the default parameters defined above, load resting state data , simulate the additon of tACS artifacts, then run artifact removal and compute the percentage errors across the spectrum for visualization.
A       = artSim(defaults);
arModes = {rbar,fbar,fastr};

pctError = table();
if ~isempty(getenv("DJ_HOST")) %&& false %[output:group:7afde059]
    %% kLab only: retrieve the data from the datajoint server,
    allRestData = ns.Epoch & 'dimension="rest"' & (ns.Experiment & 'paradigm ="rvsRestingState"');    
    sCntr =0;
    global ft_default
    ft_default.showcallinfo = 'no';
    for thisRestData = fetch(allRestData)'
        sCntr = sCntr+1;
        % Extract a Fieldtrip dataset
        ft = ephys.fieldtrip.dataset(thisRestData,dimension="rest");
        % Add the data to the ArtSim object (focus on electrodes near Oz)
        eeg = addData(A,ft,format="FIELDTRIP",channel= [125 126 138 139 117]);          %[output:3d7e592a]
        % Save the data for future use without database access.
        save(thisRestData.subject  + ".rest.mat","eeg");
        for a = 1:numel(arModes)              
            thisPctError = cleanAndEstimatePower(A,arModes{a},ar,thisRestData.subject);
            pctError = [pctError;thisPctError]; %#ok<AGROW>
        end
    end
else
    % Read subject.rest.mat data from ../data created by the code abobe
    subjects = dir('*.rest.mat');
    for sCntr = 1:numel(subjects)
        load(fullfile(subjects(sCntr).folder,subjects(sCntr).name),"eeg");
        % Add the data to the ArtSim object (focus on electrodes near Oz)
        readData(A,eeg,format="RAW",channel= [125 126 138 139 117]); 
        for a = 1:numel(arModes)               
            thisPctError = cleanAndEstimatePower(A,arModes{a},ar,extractBefore(subjects(sCntr).name,".rest.mat"));
            pctError = [pctError;thisPctError]; %#ok<AGROW>
        end
    end
end %[output:group:7afde059]
%%
%[text] ## Summarize the results
figure
tiledlayout(3,2);
for a =1:numel(arModes)
    for e = ["e5" "e10"]
        nexttile
            thisE = pctError{pctError.arMode==arMode,e};
            histogram(thisE,BinMethod="fd");
            title (arModes{a}.name +"/" + e + ": " + num2str(median(thisE),2) + "% SD: " + num2str(std(thisE),2) + "%");
    end
end
%%
%[text] ## Helper Functions
%[text] Simple example analysis; remove artifacts, then estimate the power at 10 Hz and 5 Hz, return the errors as a percentage.
function pctError = cleanAndEstimatePower(o,step,ar,subject,frequencyResolution)
arguments
    o (1,1) artSim
    step struct
    ar struct
    subject string
    frequencyResolution double =0.1
end 
o.sim;
arMode =step.name;
ar=namedargs2cell(ar);
arResult = artifactRemoval(o.vContaminated,step,ar{:},vRecordTcs =o.vRecordTcs);
r = analyzeEEG(o,arResult,freqRes=frequencyResolution,axs=[]);
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
%[output:3d7e592a]
%   data: {"dataType":"text","outputData":{"text":"Using 240.0 s of FIELDTRIP data.\nUsing 240.0 s of FIELDTRIP data.\nUsing 240.0 s of FIELDTRIP data.\nUsing 240.0 s of FIELDTRIP data.\nUsing 240.0 s of FIELDTRIP data.\nUsing 240.0 s of FIELDTRIP data.\n","truncated":false}}
%---
