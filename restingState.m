%[text] %[text:anchor:T_360F0929] # EEG Artifacts
%[text] This script investigates how tACS stimulation artifacts affect EEG/LFP analysis and shows results for FASTR, FBAR, and RBAR algorithms.
%[text:tableOfContents]{"heading":"Table of Contents"}

%%
%[text] %[text:anchor:H_D32E8D0B] ## Setup Simulation
%[text] First define the stimulation and recording parameters.
startup; % Make sure all toolboxes are on the path

defaults.reproducible              = true;    % Fix the random number generator to get reproducible results
%[text] Stimulation and Recording Equipment:
defaults.stimulatorSamplingRate    = 1e3;     % 1 kHz 
defaults.currentResolution         = 1e-6;    % 1 muA

defaults.recordingSamplingRate     = 5e3;     % 5kHz
defaults.adcLinearRange            = 200e-3;     % [V] Range over which the amplifier is linear
defaults.adcNonlinearity           = 100;       % linear inside the range,saturation outside
defaults.adcRange                  = 200e-3;     % The range of the ADC; bits are distributed evenly across this +/- range.
defaults.adcBitDepth               = 24;       % Number of bits in the ADC
defaults.adcGain                    = 1;        % Amplifier gain


defaults.highPass                  = 0.01;    % Analog filter of the recording amplifier
defaults.lowPass                   = 100;
%[text] Experiment:

defaults.duration                  = 4*60;   % 4 min of data
defaults.tacsAmplitude              = 1e-3;    % 1 mA 
defaults.resistance                 = 15;     % 15 mV artifact per mA current
defaults.tacsFrequency              = 5;      % 5 Hz tACS
%[text] And setup the base parameters for artifiact removal
arDefaults.slide                     = false;   % No sliding segment
arDefaults.slack                     = 5e-3;    % 5 ms to realign segments
arDefaults.nrSegsPerWindow           = 10;      % 10 segments to average artifact
arDefaults.pcaNr                     = 2;       % First two PC are removed
arDefaults.pcaVarExplained           = Inf;     % Ignore explained variance, always remove 2.
arDefaults.pcaNrSegsPerWindow        = Inf;     % Use all segments (plus the tCS signal; below) to determine basis set.
arDefaults.referenceSegment          = 1;       % Align to first segment
arDefaults.recordingSamplingRate      = defaults.recordingSamplingRate; % Match to the recording
arDefaults.segmentDuration          = 3;

%%
%[text] %[text:anchor:H_DDD99D7B] ## Recovering the neural background with FASTR
%[text] Create the artSim object using the default parameters defined above.
A                   = artSim(defaults);
ar                          = arDefaults;
ar.tacsFrequency            = A.tacsFrequency;
arModes = ["RBAR" "FBAR" "FASTR"];
haveDJ = true;
sCntr =0;
if haveDJ %[output:group:70b5a62b]
    allRestData = ns.Epoch & 'dimension="rest"' & (ns.Experiment & 'paradigm ="rvsRestingState"');
    pctError = table();
    for thisRestData = fetch(allRestData,'LIMIT 1000')'
        sCntr = sCntr+1;
        % Extract a Fieldtrip dataset
        ft = ephys.fieldtrip.dataset(thisRestData,dimension="rest");
        % Add the data to the ArtSim object (focus on electrodes near Oz)
        readData(A,ft,format="FIELDTRIP",channel= [125 126 138 139 117]);  %[output:2c917c69] %[output:91a41098] %[output:2136ce8d] %[output:61ad4bf2] %[output:1d0d7ab6] %[output:5aba44f9] %[output:2263ff40] %[output:05358d77] %[output:12910a9d] %[output:6e96ea19] %[output:574ff6b2] %[output:81ffdf67] %[output:08fdaa23] %[output:0502f6f1] %[output:0f3e534a] %[output:7eefe2ce] %[output:0b6286ef] %[output:7eea4e94] %[output:326d6b68] %[output:56d6a53d] %[output:95dab388] %[output:0b4eca9a] %[output:3437fbc3] %[output:39315419] %[output:19cc06d9] %[output:59e2592f] %[output:4df70499] %[output:0b53f928] %[output:564d0868] %[output:36c5b67b] %[output:17f8f3b4] %[output:5e791eed] %[output:95bf8e5d] %[output:105b4d40] %[output:6ed9e5b3] %[output:81b3658d] %[output:33f9cd8d] %[output:2d13d719] %[output:98de3323] %[output:1eb33065] %[output:4a3ac4f1] %[output:7c28e181] %[output:437e8855]
        eeg = A.lfp;
        save("../data/" + thisRestData.subject  + ".rest.mat","eeg");
        for arMode = arModes               
            r  = simPlotEEG(A,ar,arMode,freqRes =1,axs=[]); %[output:6eb2c434] %[output:04075bb0] %[output:8934818e] %[output:3eee7ea0] %[output:7a077834] %[output:5d1aad7e] %[output:22eda133] %[output:43ace6f6] %[output:6e9200ce] %[output:398c6bfe] %[output:4e25642a] %[output:25bac2dc] %[output:04903588] %[output:8d9ce638] %[output:6810eb80] %[output:8a981980] %[output:377ff7d8] %[output:672f4dd4] %[output:0c494ea7] %[output:92af5644] %[output:2d389d69] %[output:9c5a2498] %[output:5b3b2fa3] %[output:54650cea] %[output:674199b4] %[output:819cc543] %[output:89ba82e2] %[output:97f08a82] %[output:4684e91e] %[output:10acab80] %[output:116bb5bc] %[output:333a26af] %[output:7f0086ae] %[output:5b5554f6] %[output:5d7bfc15] %[output:4732e02d] %[output:04662dce] %[output:5251a5bf] %[output:5810b6dc] %[output:5378bc5d] %[output:54a20a44] %[output:988a7bc6] %[output:738c1ce4] %[output:1f639db1] %[output:27d97f6d] %[output:0449f221] %[output:0e2657db] %[output:6e20f8c2] %[output:6c3b7baf] %[output:9c5f9cef] %[output:967dafb4] %[output:984ec33d] %[output:271bdb0d] %[output:7fd4857a] %[output:5c0b757d] %[output:2bbdee8a] %[output:24c7b30b] %[output:7624a85d] %[output:6dc5a9e1] %[output:82e5210c] %[output:41d2d5d0] %[output:633a70fb] %[output:03f7d190] %[output:9f9cac9c] %[output:50c3557c] %[output:087b607a] %[output:446da65e] %[output:55089fb0] %[output:9ffef0b8] %[output:547ab557] %[output:012e47ed] %[output:38a66f3c] %[output:21feb8d6] %[output:3c17935f] %[output:661b339b] %[output:45ddae65] %[output:2ceee33e] %[output:948d86da] %[output:3f746cac] %[output:1a452905] %[output:9a709239] %[output:44081acc] %[output:0e7346aa] %[output:39d8b98b] %[output:8da3da6d] %[output:0fab5a00] %[output:931d4ea5] %[output:9baab387] %[output:2818be8e] %[output:2a557325] %[output:9546c36c] %[output:35d3eb22] %[output:2caa70ba] %[output:1d62cc69] %[output:02e537c4] %[output:06d28d2b] %[output:4248c4e7] %[output:3ea5b68b] %[output:554d52a6] %[output:6187699f] %[output:037c7147] %[output:61df3247] %[output:132a8ca0] %[output:69050d9b] %[output:80e182ca] %[output:11eaec73] %[output:92f61d11] %[output:6eaef428] %[output:27239632] %[output:6355fed7] %[output:4ba22f5a] %[output:63daa423] %[output:50c9e20c] %[output:40c9de46] %[output:9a68a6b6] %[output:9d798a2b] %[output:191d7ef2] %[output:3bcaf8d0] %[output:5964be11] %[output:32cb6b00] %[output:0ca1a68a] %[output:4481ff8f] %[output:74331512] %[output:5145b9bd] %[output:4bf2035e] %[output:2ae281d3] %[output:97da9616] %[output:72857080] %[output:60a9d174]
            s10 = r.amplitude(r.frequency>9 & r.frequency <11,1:2);
            s5 = r.amplitude(r.frequency>4 & r.frequency <6,1:2);            
            thisPct5 = mean(100*abs((s5(:,2)-s5(:,1))./s5(:,1)),"omitmissing");
            thisPct10 = mean(100*abs((s10(:,2)-s10(:,1))./s10(:,1)),"omitmissing");        
            pctError = [pctError; table(string(thisRestData.subject),arMode,thisPct5,thisPct10,'VariableNames',{'subject','arMode','e5','e10'})];            %#ok<AGROW>
        end
    end
else
    % Read subject.rest.mat data from ../data

end %[output:group:70b5a62b]
%%
figure %[output:70ab0732]
tiledlayout(3,2); %[output:70ab0732]
for arMode =arModes %[output:group:5f925c9f]
    for e = ["e5" "e10"]
        nexttile %[output:70ab0732]
            thisE = pctError{pctError.arMode==arMode,e};
            histogram(thisE,BinMethod="fd");
            title (arMode +"/" + e + ": " + num2str(median(thisE),1) + "% SD: " + num2str(std(thisE),1) + "%");
    end
end %[output:group:5f925c9f]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:2c917c69]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 2 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:6eb2c434]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 3.8%)\nFBAR (r=0.86 rmse = 31\/0 muV ae = 14.3%)\nAllocating memory...\n","truncated":false}}
%---
%[output:04075bb0]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:8934818e]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 1 sec.\nFASTR (r=0.86 rmse = 18\/0 muV ae = 33.1%)\n","truncated":false}}
%---
%[output:91a41098]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 2 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:3eee7ea0]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.02 rmse = 31\/0 muV ae = 10.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:7a077834]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:5d1aad7e]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.04 rmse = 21\/0 muV ae = 31.8%)\n","truncated":false}}
%---
%[output:2136ce8d]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:22eda133]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.01 rmse = 31\/0 muV ae = 22.0%)\nAllocating memory...\n","truncated":false}}
%---
%[output:43ace6f6]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:6e9200ce]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.02 rmse = 21\/0 muV ae = 108.2%)\n","truncated":false}}
%---
%[output:61ad4bf2]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:398c6bfe]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 1.1%)\nFBAR (r=0.79 rmse = 31\/0 muV ae = 9.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:4e25642a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:25bac2dc]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.74 rmse = 15\/0 muV ae = 50.3%)\n","truncated":false}}
%---
%[output:1d0d7ab6]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 1 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:04903588]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.1%)\nFBAR (r=0.25 rmse = 30\/0 muV ae = 9.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:8d9ce638]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:6810eb80]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.27 rmse = 21\/0 muV ae = 41.4%)\n","truncated":false}}
%---
%[output:5aba44f9]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:8a981980]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.15 rmse = 31\/0 muV ae = 13.2%)\nAllocating memory...\n","truncated":false}}
%---
%[output:377ff7d8]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:672f4dd4]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.15 rmse = 21\/0 muV ae = 29.1%)\n","truncated":false}}
%---
%[output:2263ff40]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:0c494ea7]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=-0.01 rmse = 31\/0 muV ae = 18.4%)\nAllocating memory...\n","truncated":false}}
%---
%[output:92af5644]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:2d389d69]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.00 rmse = 21\/0 muV ae = 101.9%)\n","truncated":false}}
%---
%[output:05358d77]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:9c5a2498]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.02 rmse = 31\/0 muV ae = 27.1%)\nAllocating memory...\n","truncated":false}}
%---
%[output:5b3b2fa3]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:54650cea]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.03 rmse = 21\/0 muV ae = 73.3%)\n","truncated":false}}
%---
%[output:12910a9d]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:674199b4]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.6%)\nFBAR (r=0.75 rmse = 31\/0 muV ae = 11.3%)\nAllocating memory...\n","truncated":false}}
%---
%[output:819cc543]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:89ba82e2]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.70 rmse = 15\/0 muV ae = 58.6%)\n","truncated":false}}
%---
%[output:6e96ea19]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:97f08a82]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.11 rmse = 31\/0 muV ae = 18.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:4684e91e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:10acab80]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.14 rmse = 21\/0 muV ae = 50.3%)\n","truncated":false}}
%---
%[output:574ff6b2]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:116bb5bc]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.9%)\nFBAR (r=0.72 rmse = 31\/0 muV ae = 15.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:333a26af]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:7f0086ae]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.72 rmse = 20\/0 muV ae = 41.0%)\n","truncated":false}}
%---
%[output:81ffdf67]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:5b5554f6]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.21 rmse = 31\/0 muV ae = 6.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:5d7bfc15]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:4732e02d]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.21 rmse = 21\/0 muV ae = 39.9%)\n","truncated":false}}
%---
%[output:08fdaa23]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:04662dce]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=-0.03 rmse = 31\/0 muV ae = 17.2%)\nAllocating memory...\n","truncated":false}}
%---
%[output:5251a5bf]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:5810b6dc]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=-0.02 rmse = 21\/0 muV ae = 78.8%)\n","truncated":false}}
%---
%[output:0502f6f1]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:5378bc5d]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.06 rmse = 31\/0 muV ae = 13.9%)\nAllocating memory...\n","truncated":false}}
%---
%[output:54a20a44]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:988a7bc6]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.06 rmse = 21\/0 muV ae = 47.6%)\n","truncated":false}}
%---
%[output:0f3e534a]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:738c1ce4]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.08 rmse = 30\/0 muV ae = 12.3%)\nAllocating memory...\n","truncated":false}}
%---
%[output:1f639db1]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:27d97f6d]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.09 rmse = 21\/0 muV ae = 32.3%)\n","truncated":false}}
%---
%[output:7eefe2ce]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:0449f221]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.10 rmse = 31\/0 muV ae = 14.0%)\nAllocating memory...\n","truncated":false}}
%---
%[output:0e2657db]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:6e20f8c2]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.09 rmse = 21\/0 muV ae = 48.7%)\n","truncated":false}}
%---
%[output:0b6286ef]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:6c3b7baf]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.1%)\nFBAR (r=0.18 rmse = 31\/0 muV ae = 12.9%)\nAllocating memory...\n","truncated":false}}
%---
%[output:9c5f9cef]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:967dafb4]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.17 rmse = 21\/0 muV ae = 56.5%)\n","truncated":false}}
%---
%[output:7eea4e94]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:984ec33d]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.00 rmse = 31\/0 muV ae = 22.0%)\nAllocating memory...\n","truncated":false}}
%---
%[output:271bdb0d]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:7fd4857a]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.02 rmse = 21\/0 muV ae = 80.2%)\n","truncated":false}}
%---
%[output:326d6b68]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:5c0b757d]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 5.1%)\nFBAR (r=0.96 rmse = 30\/0 muV ae = 12.3%)\nAllocating memory...\n","truncated":false}}
%---
%[output:2bbdee8a]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:24c7b30b]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.93 rmse = 27\/0 muV ae = 63.5%)\n","truncated":false}}
%---
%[output:56d6a53d]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:7624a85d]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.10 rmse = 31\/0 muV ae = 18.6%)\nAllocating memory...\n","truncated":false}}
%---
%[output:6dc5a9e1]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:82e5210c]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.09 rmse = 21\/0 muV ae = 53.1%)\n","truncated":false}}
%---
%[output:95dab388]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:41d2d5d0]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.02 rmse = 31\/0 muV ae = 17.4%)\nAllocating memory...\n","truncated":false}}
%---
%[output:633a70fb]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:03f7d190]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.03 rmse = 21\/0 muV ae = 74.8%)\n","truncated":false}}
%---
%[output:0b4eca9a]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:9f9cac9c]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.02 rmse = 31\/0 muV ae = 21.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:50c3557c]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:087b607a]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.03 rmse = 21\/0 muV ae = 106.2%)\n","truncated":false}}
%---
%[output:3437fbc3]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:446da65e]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.10 rmse = 31\/0 muV ae = 12.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:55089fb0]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:9ffef0b8]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.10 rmse = 21\/0 muV ae = 33.1%)\n","truncated":false}}
%---
%[output:39315419]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:547ab557]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.12 rmse = 31\/0 muV ae = 8.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:012e47ed]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:38a66f3c]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.11 rmse = 21\/0 muV ae = 32.6%)\n","truncated":false}}
%---
%[output:19cc06d9]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:21feb8d6]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.09 rmse = 31\/0 muV ae = 15.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:3c17935f]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:661b339b]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.08 rmse = 21\/0 muV ae = 64.5%)\n","truncated":false}}
%---
%[output:59e2592f]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:45ddae65]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.09 rmse = 30\/0 muV ae = 15.4%)\nAllocating memory...\n","truncated":false}}
%---
%[output:2ceee33e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:948d86da]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.08 rmse = 21\/0 muV ae = 48.1%)\n","truncated":false}}
%---
%[output:4df70499]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:3f746cac]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.1%)\nFBAR (r=0.13 rmse = 31\/0 muV ae = 15.4%)\nAllocating memory...\n","truncated":false}}
%---
%[output:1a452905]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:9a709239]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.12 rmse = 21\/0 muV ae = 49.4%)\n","truncated":false}}
%---
%[output:0b53f928]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:44081acc]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.21 rmse = 31\/0 muV ae = 15.1%)\nAllocating memory...\n","truncated":false}}
%---
%[output:0e7346aa]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:39d8b98b]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.20 rmse = 21\/0 muV ae = 46.5%)\n","truncated":false}}
%---
%[output:564d0868]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:8da3da6d]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.05 rmse = 31\/0 muV ae = 18.5%)\nAllocating memory...\n","truncated":false}}
%---
%[output:0fab5a00]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:931d4ea5]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.06 rmse = 21\/0 muV ae = 51.3%)\n","truncated":false}}
%---
%[output:36c5b67b]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:9baab387]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.03 rmse = 30\/0 muV ae = 25.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:2818be8e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:2a557325]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.03 rmse = 21\/0 muV ae = 60.0%)\n","truncated":false}}
%---
%[output:17f8f3b4]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:9546c36c]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.34 rmse = 30\/0 muV ae = 5.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:35d3eb22]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:2caa70ba]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.37 rmse = 21\/0 muV ae = 10.9%)\n","truncated":false}}
%---
%[output:5e791eed]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 3 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:1d62cc69]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.12 rmse = 31\/0 muV ae = 14.9%)\nAllocating memory...\n","truncated":false}}
%---
%[output:02e537c4]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:06d28d2b]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.11 rmse = 21\/0 muV ae = 41.1%)\n","truncated":false}}
%---
%[output:95bf8e5d]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:4248c4e7]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.04 rmse = 31\/0 muV ae = 23.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:3ea5b68b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:554d52a6]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.05 rmse = 21\/0 muV ae = 74.9%)\n","truncated":false}}
%---
%[output:105b4d40]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:6187699f]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.3%)\nFBAR (r=0.17 rmse = 30\/0 muV ae = 11.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:037c7147]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:61df3247]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.17 rmse = 21\/0 muV ae = 28.0%)\n","truncated":false}}
%---
%[output:6ed9e5b3]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:132a8ca0]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.03 rmse = 30\/0 muV ae = 21.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:69050d9b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:80e182ca]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.03 rmse = 21\/0 muV ae = 63.1%)\n","truncated":false}}
%---
%[output:81b3658d]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:11eaec73]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.03 rmse = 30\/0 muV ae = 8.9%)\nAllocating memory...\n","truncated":false}}
%---
%[output:92f61d11]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:6eaef428]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.05 rmse = 21\/0 muV ae = 38.3%)\n","truncated":false}}
%---
%[output:33f9cd8d]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:27239632]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 7.5%)\nFBAR (r=0.96 rmse = 31\/0 muV ae = 14.9%)\nAllocating memory...\n","truncated":false}}
%---
%[output:6355fed7]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:4ba22f5a]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.92 rmse = 63\/0 muV ae = 85.7%)\n","truncated":false}}
%---
%[output:2d13d719]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:63daa423]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.10 rmse = 31\/0 muV ae = 12.3%)\nAllocating memory...\n","truncated":false}}
%---
%[output:50c9e20c]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:40c9de46]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.09 rmse = 21\/0 muV ae = 34.9%)\n","truncated":false}}
%---
%[output:98de3323]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:9a68a6b6]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 1\/0 muV ae = 35.4%)\nFBAR (r=1.00 rmse = 31\/0 muV ae = 39.6%)\nAllocating memory...\n","truncated":false}}
%---
%[output:9d798a2b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:191d7ef2]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.97 rmse = 147\/0 muV ae = 187.3%)\n","truncated":false}}
%---
%[output:1eb33065]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:3bcaf8d0]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.2%)\nFBAR (r=0.14 rmse = 31\/0 muV ae = 8.1%)\nAllocating memory...\n","truncated":false}}
%---
%[output:5964be11]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:32cb6b00]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.13 rmse = 21\/0 muV ae = 38.8%)\n","truncated":false}}
%---
%[output:4a3ac4f1]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:0ca1a68a]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.1%)\nFBAR (r=0.13 rmse = 31\/0 muV ae = 8.9%)\nAllocating memory...\n","truncated":false}}
%---
%[output:4481ff8f]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:74331512]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.13 rmse = 21\/0 muV ae = 50.3%)\n","truncated":false}}
%---
%[output:7c28e181]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:5145b9bd]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 0.1%)\nFBAR (r=0.16 rmse = 31\/0 muV ae = 17.8%)\nAllocating memory...\n","truncated":false}}
%---
%[output:4bf2035e]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:2ae281d3]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.14 rmse = 21\/0 muV ae = 58.3%)\n","truncated":false}}
%---
%[output:437e8855]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds and required the additional allocation of an estimated 0 MB\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:97da9616]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 1\/0 muV ae = 11.2%)\nFBAR (r=0.99 rmse = 31\/0 muV ae = 31.7%)\nAllocating memory...\n","truncated":false}}
%---
%[output:72857080]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0.0482"}}
%---
%[output:60a9d174]
%   data: {"dataType":"text","outputData":{"text":"\nStage 1 of 2: Slice Alignment\n....25% ....50% ....75% ....100%\n\nStage 2 of 2: Artifact Subtraction\n25% 50% 75% 100%\nFASTR Finished in 0 min 0 sec.\nFASTR (r=0.97 rmse = 62\/0 muV ae = 195.6%)\n","truncated":false}}
%---
%[output:70ab0732]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tvQu4lUX5\/v+IIqaIKB7wUGEmiGICJWhqHlILPKap5QHD8pxCoaKYIiiiFSKplSKEqagoKBEJpaZ5yBQl0lKQsDLUnxYmSYYc\/P8\/83V2ay\/XXjN7HfZ7WPdcFxew17wzz9zPYe55nnnXXuf9999\/39SEgBAQAkJACAgBIZAhBNYRgcmQtiSqEBACQkAICAEh4BAQgZEhCAEhIASEgBAQAplDQAQmcyqTwEJACAgBISAEhIAIjGxACAgBISAEhIAQyBwCIjCZU5kEFgJCQAgIASEgBERgZANCQAgIASEgBIRA5hAQgcmcyiSwEBACQkAICAEhIAIjGxACQkAICAEhIAQyh4AITOZUJoGFgBAQAkJACAgBERjZgBAQAkJACAgBIZA5BERgMqcyCSwEhIAQEAJCQAiIwMgGhIAQEAJCQAgIgcwhIAJTocr+9Kc\/2aGHHtr09Prrr28f+9jH7LTTTrMvf\/nL7ufFfTp06GDbbLONHX\/88fb1r3+92czXXXedjR8\/3v3s\/PPPtzPPPLPp88mTJ9sVV1zR9P+OHTtaz5497cILL7Q+ffo0G+eb3\/ymbbLJJjZmzJiyK3vsscds0KBBH+rzq1\/9ynbYYYcoVG677Tb7yU9+YkuXLjXWz3NDhw61fffd1z1\/55132jXXXGP\/\/e9\/7bDDDnNrWGeddWz58uV24IEH2pVXXun+bqktXLjQ9fn9739vq1atsq233tr69u1r3\/ve99wjhbgw\/8Ybb2w9evRw+A4cODBqDb7T7Nmz7YYbbrAlS5bY5ptv7vQzePDgVo0R0\/n++++3cePG2SuvvOJs4eyzz26yl9dee80uuugimz9\/vv373\/+27bff3h588MGmYavFM0Y+9WkZgd12283phbbeeuvZlltu6WIA\/rruuuu6nxf3wZb2339\/+853vmMf+chHmgb\/3e9+Z1\/96lfd\/\/fZZx+75ZZbmj6LiS2+s\/fBQjtpaQX4y7Rp02zRokWuy4QJE5xf+vb222\/bqFGjnM2tWbPGPve5z7n\/b7HFFtFm8cILL7hnsOFOnTrZl770JbvgggscXqUav0sYOe644w5jfuLapZde2hTXfvSjH9msWbPs73\/\/u4sdO++8s4t74Jy1GNNabMqtHdywKbDZYIMN7Nvf\/rZ95StfcZiA\/cknn2wPPPCAs9E8NxGYCrVbGGQIUGw+L774ohttzpw51r1792YEhmD15ptvGgQBR5w6dartscceTbMfdNBB9uc\/\/9n9f6eddrJf\/OIXTZ\/5jZqAAGEhADFfly5d7Omnn27qR9D59Kc\/bd\/\/\/vfLEgMe8ASGTXK77bZrGmPs2LFuY41pI0eOdDLjJOCBXBCJp556ytq3b2+f+tSn3Abdr18\/O\/HEEx3Zgdxcfvnljijw\/5ba2rVrXQB99dVXHTHaaqut7G9\/+5v94x\/\/cHPRPC4QOvpACpYtW+Y+Y95hw4bFLMOeeOIJO+GEExxxIQigI9Zy\/fXXt5oIvfPOO3bEEUfYVVddZbvvvnuz+b3NgNF+++3n5qU\/m9BnP\/tZI8AhN\/jzWSGBeffdd6vCMwoIdSqLgCcn\/M1hBDunYc\/YD833+fznP+8INZsIOj733HMdufeNTRq90yA\/v\/3tb5390WJiix\/n1FNPdQenSy65JKg9Njn8CT96\/fXXP0RgTj\/9dGf7xB82RQ4O2PBdd90VHJsOK1eudD5LnIOU4eMcbsr5IjEA\/CBJxAvIE7j95je\/cQcxDiKQRnD9wx\/+4IjMpptuao888ojDLSsxphJsyq3dkxRsiD2AWMXPNtxwQ0caDznkEHeYznsTgalQw4VBBkeFEeO8OCzZFAyouA9TceL54x\/\/6LIRZApobFz0h5BAbtik586dazvuuGOzjRoiwCmcOQgQNE5y\/oT0zDPPuDH5m02dTRAy89JLLzmnZ3M+44wzrF27dk0EpjjbUwgHpx2yJ4VztAQXWZXevXu7jz3z33XXXe2yyy6zvffe2yBoN910k9uUwYCMB\/9uqfk1EtSZH1xoECafIfIExuMC6eE0B\/40giFzkE359a9\/\/aGA7eeGiN16661Nmwzy4\/zIP3PmTNetHJaFa\/A4eLJW+BkbyH333Wff+ta37JxzzrGf\/vSnDh\/sZsqUKU1dCURHH310MwLDJlgNnhWauR4rQMCTE3\/4OO+882zGjBnO53yGtLgPGcSbb77ZZWp+8IMfuNFWr17tDi+QbcgCB5\/Ro0c7kk+LiS30e++991xGkpM68QBighzPPvus83EOCxAbDj6FjdP5o48+2swfXn75ZTvggANso402chsiB5C99trLkRHWiG+H\/AiiQwaRtYERh4AvfvGLLhbNmzfPHW6KG8QdMkUGpn\/\/\/kYGmcMbWRtiFb6Pf\/u4SB9i7fTp0118rMYn2jLGVIJNubW\/9dZb9o1vfMPI6IIt2RgIHnF14sSJ7hCNDvPeRGAq1HBxkMEZICFsYKEMDA6NoXFyol199dV244032rHHHuv+T5q38NRSuFHj6GzGGC+E58knn2xKX1OaYPODlUOKDj\/8cHeSOumkk1xQgtiQfmVz9hkYAhblmW233dadIk855ZQmRGIIDOM8\/vjjbl5OpJw8cSDaD3\/4w6agzSZNiYYguMsuu9jw4cPLIr9ixQr7zGc+4051BEECF8HOkzoeLiYw\/IzgBmEiQwUxIViHAq\/fZChnQYCQG1kJuGwuISyZl5MhjdMidvDd7363KcMGweTEPmDAAKMsxobGZuGJCp8TrHwrRWCqxbNCM9djLRCYXr16OR\/E5iGhvhxbnIGBRBMTOCFjFzSyB5Qnu3Xr5jZsiFBhpiMmtjAOvkcGhkwJpOgLX\/iCy7AQRyA3kGVITHGmsxSBYSMk5pDR4Dma74d\/cPgJ+ZE\/CIDLiBEj3BgQNGTxh4lCg2ITJmNM41BHiQ0yRomYQw6+WNjwaWIA5IyMFZnfrMSY1mJT7HjFa99ss83srLPOcpkqGnaEvogrHFp9GT\/vDiwCU6GGi+vUDMOGR\/33uOOOc6OW6sPPCXJs8j5lzIZL4GFjI4sAwfj4xz\/uiAqt+A4MP+vatavbZAvvwOD0kBaCGoycUxBGDinhlMMpkHskEA6CH4EXMvGvf\/3LOQKnHIz\/qKOOcvMSyAiMBN7C+n0hZJAm5KCRvuQE5tPp\/IznIUg8\/8tf\/tKdCDnRUWLhxAAWnBqLT4k8S1Al1f7Pf\/6zaUpOhWRyGK8UgaEj+EPYPAlkrawfEvTRj370Qxrn5ApubDQ0SJDP+LA5QEbKYckaKRm21PyJnZMk2TX\/f8gMpIZUOFky31oiMNXiWaGp67EPECi83+JBoTSMz\/s7HqX6QFA5oPgMpc\/csNFjoxB1bAg75SARE1uYH2JBRnLSpEnuQERWDxtjU6dh08SVhx9+uOmwxM9LEZjbb7\/d+abPZtKP2MFhjMwhG2TIj5gfOQrLZcQn7raQgSguqSI7mVn8jQwQjbtA4AlR8SU2fs6BhjiGX3t5vA6yEGNai02h05VbO6VlbI9sC1k8yujgA4ZkzyhDxZbSs+joIjAVaq0wyEAcKDGw0WIwnLZoxScpjI3UKE5O4OOiLVmRY445xm3+\/BsCQ1qYzMO9997rNni\/UUN4CIKkf\/m88OSHsRJ8KD2xmX7ta19rYufFS+S0Q0aAk4xv\/oRAQCHYtqZBUNjoOVX+5z\/\/aUo5F46BvGQ4cKYFCxa4WjskjiAOQfIntuJ5weO5555zpOvHP\/6xu0\/gT4QtZWAgORAWn4GJWQvEgtQ1JIbME3IRFMi+EDj9SacUlmS5WA+N9RNACFiQQxobFCcmn4Fhw+HeVGszMLXAMwYL9SmNgCcnlGs4ES9evNhlUSgzcm+DVlhCYsOmTMgdD3+fCT\/g59ix39TJbpDJgfxz+IiJLcx18MEHu7IT2R98A6JdqhWXM8tlYJCfuEMrzsCE7MLHEA5gZHppXKonPlSTgcE38UcOPNzTaSl7m+YY01psPNaxa6dcx8GTGEbpjYMa+sM+OPCVe1kipNc0fy4CU6F2iskJJx3KJziRDxil7sBAbngzx58wvGGXEoOU7cUXX\/yhTAOnIk5HkB42VjZHgg7ZEzZ6Gs9RbsKACXS+kc6lnANTJ8vjSYyXg03WZ1S44U7wgWBAeIobJysu2tG4QMzFVEpprO\/II49s1p17KRAvymOQN54j6BKwwYxAX9gI8Jz4cDx\/uvUnQl8GKyYwyMA9g+I7MJTZ0A8ZGE64xQ2SRMbFZ11weDJEZMaQK4Rl4RsW5e7AcPeFzc6fIP0dmOK3UMplYLzsrcWzQjPXYwUIFJMTMhwQ3MJNtfgOjD+gkGVjk8F3ybqUahDnn\/\/85x86+JSKLf7+hs+u4KtDhgxxJRn+LvR3ykI+21tITArfQip1B2bPPfd0hzLum5BJCfkR9\/M4iPjY5u\/AQO4gaMQa5MRfwI7GfZk33nijicyBDZlXfzeP0iwEDfnIKhe\/vVmIY2t9oi1jTCXYtGbtxHl0RGz55Cc\/6TLXEBjuCGGfZNDy2ERgKtRqKXJC2o4UqE\/DFvZh0+ZCLGUUUoKwZLIRODBBglOZJwMQA1Kl1HjJ7HDBkzJLYXqXTA\/3M3y6lr+5W0N2gvb888+7t2EIHpzwyCZQsiBDQGoWEgBBwOhxZOryNH8BmX+H7sAQrMkIkSInkPPHZyMKiQLBlnssnDgZkyBHMEQOZCMbAXaFjYwShIPSEvdeSBOTuSHdDFljXk9guMfziU98wv761782lYEK7xCFavfMhRwEbC5WEkDBC3nBJ4RlodzlCAzjkK0DI8guhI47MxAZyBJycB8Ke0Af6BPyCdEko+NbJXhWaOZ6rAUCg9\/6gwR64uCArxXegYE04GPoC3uldOo3aC6iF15i928nkqkgU+u\/ooEXBGjFsYUyJJk8\/\/o0PoyPQXYoofImG\/7A58QSSAPPcMEXu8PWiDmc1PEPXl8mc8IFdjKH2CjkC0J09913OxlCfoTcHI6wX+5iUBbloOR90fuGv1vGmN6HKYmz2RIfwRH758UDDkJkXvgZn\/vGYcDfn+FnlfhEW8aYSrCJXTvxiv0B3VFaJ7ZAWjkccvgsvH+VN4cWgalQo6UIDE7EZsxmy6mFrEXhd8VwCoOUsClBIAgsnC4gLmzoPsuBsRMg2dwIOsxVTGD8iYtnCUic4nkFmiDmG8GT75bx9ysImFzw444Kxk4AhARBrPiMNK2\/\/xJDYNhUOVlxGQ85ICcw\/eJad\/F301De4ZVSCAmBnaxJ8XdNkJXhhMjaCMRkgpCRIOuzO4V3gyAckB1S1pBFf2EyJvCCM3oAJ+bB+T25jMEylsDQjxSv\/x4YgjbB3d+Z4pRJ8CluhWl9PqsEzwrNXI+VITCQaUqukAyfMSi8A8PnXLSnhEgWjw2Zf5MFJRNXuCH7cg0HEeJDMYEpji1cduWQQObUN+yHzCEkBZvmc4gxpWayiz4DWKxUnzHGj+nL3TtiGDEF4oSdxvgRfTjEMAa+DbHjEEWswz9LERj\/5qD\/HhhiCJkWT07AglJdcSN7W5hZrsQn2jLGVIJNzNpZAzGDspr\/Th9eCOCQSGYLO+KqQuF1gTw5dcMRGJycMgns1F9uwnHZUHGUzp07u02\/3KXMtBkAJQc2QQIXQUNNCAiB0gjkwf8h2dyT42TdKG+byJ6FQCkEGorAcCLicivpSE7znsBwyuCEAGvlgi0nezIoWWkQGC7m+u+RyIrcklMItCUCefF\/Sh+8NcTdhpbeDmxLXDWXEEgKgYYiMNy5INtCCpfmCQx3S7hcyhsF\/i0g6ts+dZqUcjSvEBACtUNA\/l87LDWSEEgDAg1DYKghU1\/mNrj\/RkwIjP8ODy64+m8u5IIsdcPi3zOUBoVJBiEgBFqPgPy\/9ZjpCSGQdgQahsBwQZVfkMjlMC5R+gwM9WQufhYSGO7IcCHWf\/FU2pUo+YSAECiPgPxfFiIE8odAwxAYbv\/zaiDN\/0ZZLr5y7yVUQuJ3nfCWkJoQyCoC\/vfTZFX+auWW\/1eLoJ7PMgJ59f+GITCFxleYgeHnvIrIK8z+Ei+vF\/OdDb7xHSP++xhKGTGv+xZ+eVSWDb0R18ea867DkA3nzWbLrUf+n7y28+5v9UC4Gszy6v8iMGbuYi8ExL9GTYDji51iCUw9jFVjCoFaIpDXAFYJRsUERv5fCYp6JksI5NX\/G5LAtNbwQsrvd9CR9o8\/\/6HFYfOavmstjuqfHAIhG05OsvTPLOzSryNJWB6BvNqwCEyE5YeUz+dHXH1\/iyPNHD6gbAkqQgR1EQJVIRCy4aoGz\/nDIezK3ZHT4SXnxpGR5YVsOCPL+JCYIjARmgspXwQmAkR1SRSBkA0nKlzKJw9hV87\/dXhJuXIbRLyQDWcVBhGYCM2FlC8CEwGiuiSKQMiGExUu5ZOHsBOBSbkCJZ77ZbflXkTJKkQiMBGaCylfBCYCRHVJFIGQDScqXMonD2EnApNyBUo8EZhGtoFqAhi4KY3cyNaTjrWHbDgdUqZTihB2IjDp1Juk+h8CIRvOKlbKwERoLqR8ZWAiQFSXRBEI2XCiwqV88hB2IjApV6DEUwamkW2gmgCmDEwjW0561h6y4fRImj5JQtiJwKRPZ5KoOQIhG84qXsrARGgupHxlYCJAVJdEEQjZcKLCpXzyEHYiMClXoMRTBqaRbaCaAKYMTCNbTnrWHrLh9EiaPklC2InApE9nkkgZGNnABwhUE8BEYGRGaUAgZMNpkDGtMoSwE4FJq+Ykl0cgZMNZRUolpAjNhZSvElIEiOqSKAIhG05UuJRPHsJOBCblCpR4KiE1sg1UE8CUgWlky0nP2kM2nB5J0ydJCDsRmPTpTBKphCQbUAlJNpATBEKbcE6WWZdlhLATgakL7Bq0hgiEbLiGU7XpUCohRcAdUr5KSBEgqkuiCIRsOFHhUj55CDsRmJQrUOKphNTINlBNAFMJqZEtJz1rD9lweiRNnyQh7ERg0qczSaQSkmxAJSTZQE4QCG3COVlmXZYRwk4Epi6wa9AaIhCy4RpO1aZDqYQUAXdI+SohRYCoLokiELLhRIVL+eQh7ERgUq5AiacSUiPbQDUBTCWkRrac9Kw9ZMPpkTR9koSwE4FJn84kkUpIsgGVkGQDOUEgtAnnZJl1WUYIOxGYusCuQWuIQMiGazhVmw6lElIE3CHlq4QUAaK6JIpAyIYTFS7lk4ewE4FJuQIlnkpIjWwD1QQwlZAa2XLSs\/aQDadH0vRJEsJOBCZ9OpNEKiHJBlRCkg3kBIHQJpyTZdZlGSHsRGDqArsGrSECIRuu4VRtOpRKSBFwh5SvElIEiOqSKAIhG05UuJRPHsJOBCblCpR4KiFl3QbuvPNOmzRpki1fvty6dOliI0eOtP79+7tlvfXWWzZ06FBbvHixde7c2caPH2\/du3dvWnI1AUwlpKxbTj7kD9lwPlbZ8iqS8v+ZwwfYkiVL8g6v1pdyBPLq\/w2TgXnooYds1113tS222MIeffRRR1iefvppa9eunY0aNcrat29vI0aMsNmzZ9vkyZNt+vTpIjApd0qJF49AXgNYLAJJ+b8ITKyG1K+eCOTV\/xuGwBQax3vvvWe9evWyZ5991jp27Gj9+vWzadOmWbdu3Wzt2rXWt29fmzNnjnXt2tU9FlK+Skj1dD2NXQsEQjZcizmyMkZb+r8ITFasIt9y5tX\/G5LA3H777TZr1iwjrbx69WpXLlq4cKHLwtAGDhxoY8aMsT59+ojA5NuvG2Z1eQ1glSiwLf1fBKYSDemZWiOQV\/9vOALzxBNP2MUXX2y33XabbbvttrZq1Srr0aNHMwIzYMAAGzt2rPXu3VsEptaepPESQSCvAay1YLa1\/4vAtFZD6l8PBPLq\/w1FYLjzMmzYMJsyZYorC\/kWU0IaMmSI8adUUwmpHi6nMWuFwIQJE4w\/jX6ZNAn\/F4GplRVrnEoRyLP\/NwyBWbBggZ199tk2ceJE69mzZzNb4I2kDh06NF3i5W2lGTNmNPUJsVcRmEpdS8+1FQIhG24rOZKaJyn\/F4FJSuOatxCBvPp\/wxCYQYMG2bx586xTp05NeuUODBd3eY2a7Ip\/jXrcuHHNSE5I+SIwChZpRyBkw2mXv1r5kvJ\/EZhqNafna4FAXv2\/YQhMNUYQUr4ITDXo6tm2QCBkw20hQ1bnCGGnL7LLqmYbR+6QDWcVCRGYCM2FlC8CEwGiuiSKQMiGExUu5ZOHsBOBSbkCJV7wq0CyCpEITITmqglgDK80cgTI6lJXBEI2XNfJMz54CDsRmIwruAHED9lwViEQgYnQXEj5ysBEgKguiSIQsuFEhUv55CHsRGBSrkCJpwxMI9tANQFMGZhGtpz0rD1kw+mRNH2ShLATgUmfziRRcwRCNpxVvJSBidBcSPnKwESAqC6JIhCy4USFS\/nkIexEYFKuQImnDEwj20A1AUwZmEa2nPSsPWTD6ZE0fZKEsBOBSZ\/OJJEyMLKBDxCoJoCJwMiM0oBAyIbTIGNaZQhhJwKTVs1JLo9AyIazipRKSBGaCylfJaQIENUlUQRCNpyocCmfPISdCEzKFSjxVEJqZBuoJoApA9PIlpOetYdsOD2Spk+SEHYiMOnTmSRSCUk2oBKSbCAnCIQ24Zwssy7LCGEnAlMX2DVoDREI2XANp2rToVRCioA7pHyVkCJAVJdEEQjZcKLCpXzyEHYiMClXoMRTCamRbaCaAKYSUiNbTnrWHrLh9EiaPklC2InApE9nkkglJNmASkiygZwgENqEc7LMuiwjhJ0ITF1g16A1RCBkwzWcqk2HUgkpAu6Q8lVCigBRXRJFIGTDiQqX8slD2InApFyBEk8lpEa2gWoCmEpIjWw56Vl7yIbTI2n6JAlhJwKTPp1JIpWQZAMqIckGcoJAaBPOyTLrsowQdiIwdYFdg9YQgZAN13CqNh1KJaQIuEPKVwkpAkR1SRSBkA0nKlzKJw9hJwKTcgVKPJWQGtkGqglgKiE1suWkZ+0hG06PpOmTJISdCEz6dCaJVEKSDaiEJBvICQKhTTgny6zLMkLYicDUBXYNWkMEQjZcw6nadCiVkCLgDilfJaQIENUlUQRCNpyocCmfPISdCEzKFSjxVEJqZBuoJoCphNTIlpOetYdsOD2Spk+SEHYiMOnTmSRSCUk2oBKSbCAnCIQ24Zwssy7LCGEXIjDlhNpjjz1s6tSpdZFbgwoBj0DIhrOKlEpIEZoLKV8lpAgQ1SVRBEI2nKhwKZ88hF2IwBxx9f0trnDm8AG2ZMmSlCMg8bKOQMiGs7o+ERgze+utt2zo0KG2ePFi69y5s40fP966d+\/epNOQ8kVgsmr+jSN3yIYbB4kPr7Se\/g9BEYFpZOtKx9rz6v8iMGY2atQoa9++vY0YMcJmz55tkydPtunTp7cpgTn++OPtySefbNHalWpORyDIqhR5DWC10Ec9\/V8EphYa0hjVIpBX\/xeBMbN+\/frZtGnTrFu3brZ27Vrr27evzZkzx7p27ersJqT8WmRgajFGtUau5\/OLQMiG87vy8Mrq6f8iMGH81aP+COTV\/xuewKxevdqVixYuXOiyMLSBAwfamDFjrE+fPqkiMCEzr0WWplEzQRMmTLAhQ4aEIE7s82r1ktcAVq1C6u3\/MQSm3Bpq4dPVYlSP59Pub\/VYc7VjlosBITvJq\/83PIFZtWqV9ejRoxmBGTBggI0dO9Z69+7tbC60eVRrmHpeCNQbgVCAq\/f8aR1f\/p9WzUiuWiKQV\/9veAKDkYRSyLU0JI0lBIRAuhCQ\/6dLH5JGCMQiIAJjZiNHjrQOHTo0XeKdNGmSzZgxIxZD9RMCQiDDCMj\/M6w8id7QCIjAfPAaNfcf\/GvU48aNs549eza0YWjxQqBREOA1avl\/o2hb68wTAiIwVWiTLM0111xj77\/\/vh1xxBF2wQUXVDFaMo+++eabLvP03HPP2brrrmv77befXXbZZU0XmltaY+i7M5JZTXjWiy66yO666y5btGiRrbfeeu6BvK0xjIJ61AKBPPh\/LXAoNQYx8fLLL3dvcxJXLrzwQjvkkEM+1PXVV1+1fffd17p06eI+IxP+yCOP1EuszI27fPlyO+uss1x85iLuvffem7k11FNgEZgK0WXjxyHvu+8+53xHHXWUc9J99tmnwhGTeYx1PP\/88y6IrFy50k499VQ76KCD7OSTT7Zyawx9d0Yyqyk\/62OPPeYCAH88gcnbGtOIex5lyov\/10s3EBdK8fyahNdee82OPvpoe+ihh2zjjTduNiUE5qSTTrIHH3ywXqJkelxi8vz58+2dd96xG264QQSmSJsiMBWaN445b948l4Gh8eV3bIpXXXVVhSOm4zHKZzgNWZlya8zaxcd3333XvU120003Wf\/+\/ZsITJ7WmA4Lagwp8ur\/tdLeueeea3vttZcdd9xxbsjTTjvNfT3FkUceKQJTAch8yenVV18tAiMCU4H1lHiEjZ406Xnnnec+5cRxxx132C233FKbCRIYZcWKFXbooYe6V8h57a6lNXKyCn13TgLil51y9OjR7rV4giiy+wxMntaYNszzLE8e\/b+W+jr22GPtnHPOacpIk7Hdcsst7cwzz\/wQgdl\/\/\/1tq622sg033NBlY0444YRaipKLsURgSqtRGZgKzZsARhs2bJj7+\/7773d3K6ZMmVLhiMk+xvdhnHHGGbbrrru63wtFa2mNEydODH53TrKraT47Kdjrr7\/epbT9F5cVEphSeszaGtOEdyPIkjf\/r7XOjjnmGHcxeu+993ZD86bX1ltv7WJMYSPbu2zZMvfZyy+\/bIMGDbIrr7wyc6X4WuNXPJ4IjAhMTW0sTynkNWvWGClfTkGXXnppE055Ka9ARviz\/vrru7VRd99mm23cr494+OGHWywFZq1MVlMD12BlEciT\/9dD1bElpOK5+Qb0Tp06ueyN2v8QEIERgampP\/gTlPRcAAAgAElEQVRLfDNnzrTNNtvMXVIbPnx45k4O\/O4nskgbbLCBO\/mss846TTiVW2NWvzujOAOTxzXW1NA1WEkE8uL\/9VIvGWnuBfpLvLzkwCVeyMmCBQvcZV7eqlm6dKltsskm1rFjR3v99ddd+Yg3BQ888MB6iZbJcUVgRGBqbrj33HOPXXvtte4XQB5++OHuLaSstWeeecZI926++ebWrl07J\/7BBx9s3BmhtbTGrH53RjGByeMas2aDWZU3D\/5fL+y5H0gMmTt3rosrkBL\/GjX3BnfYYQd3H4a7g1dccYVRwuYVai7aF5eZ6iVjVsblrVDKbLyJtOmmm9rpp59ugwcPzor4dZVTd2DqCq8GFwJCQAgIASEgBOqBgAhMPVDVmEJACAgBISAEhEBdERCBqSu8GlwICAEhIASEgBCoBwIiMPVAVWMKASEgBISAEBACdUVABKau8GpwISAEhIAQEAJCoB4IiMDUA1WNKQSEgBAQAkJACNQVARGYusKrwYWAEBACQkAICIF6ICACUw9UNaYQEAJCQAgIASFQVwREYOoKrwYXAkJACAgBISAE6oGACEw9UNWYQkAICAEhIASEQF0REIGpK7waXAgIASEgBISAEKgHAiIw9UBVYwoBISAEhIAQEAJ1RUAEpq7wanAhIASEgBAQAkKgHgiIwNQDVY0pBISAEBACQkAI1BUBEZi6wqvBhYAQEAJCQAgIgXogIAJTD1Q1phAQAkJACAgBIVBXBBqWwEyePNmuuOKKD4F70UUX2amnnmqnnXaaPfDAA02fb7bZZtanTx8bMWKEbb\/99s2eO\/744+3JJ590P5s+fbrr51vhOO+\/\/75tscUW1q9fPxs5cqT7t28vvfSSHXzwwfbUU081+3kp7Y8bN85uuOGGZh9tvvnm7tnYxvqnTZtmixYtco9MmDDBDjvssNjHo\/u98MILNmrUKJs\/f7516tTJvvSlL9kFF1xg6623nhvjvPPOc9i9+uqr7v+\/+93vmtbPM8OHD7e\/\/e1v9pnPfMbGjx\/f9NnJJ59sPXr0cPpoqYH3j370I7vzzjvtjTfesI033th22GEH++Y3v2l77723\/elPf7JDDz3UPb7uuuvaRz7yEafbL37xi\/b1r3\/d1l9\/\/eh1vvnmm3bVVVfZo48+av\/+979tjz32sMsvv9y222676DFiOobwLKfXf\/zjH\/atb33L5s2b5+T67ne\/22Sr4DR37ly79957bZ111okRpeH6hGJG8efY2yc+8QkbOnSo7bvvvjWPGdg3seTSSy8N+u5rr71mxDZ8CvvEzh988MFmMt1\/\/\/1GbHnllVdsm222sbPPPtu+\/OUvR+sZeYgjd9xxh7399tvWs2dPJ1thPCweLGTPs2fPdrFuyZIlRozDLwcPHuyGqXd8YI7ddtvN4UXbaKONbMstt7R99tnHzjrrLPfv2NZabP785z\/b1VdfbX\/4wx\/sX\/\/6l1s7sRNbIlZVu\/ZYudPer+EJTMeOHd1G6NugQYNcMPDEY5dddrFevXq5jYlNloDBhugbGxebFQZK+9rXvuac1jc\/Trdu3eyjH\/2oPfvss7ZixQo75JBD7LrrrmvqN2nSJJsxY4bhsKHmCQzOBSmgbbLJJvaDH\/wg9GjT59\/+9rfdeiAHr7\/+elUEhjXuvvvujvgVtpUrV9rnPvc5AyOcniC0dOlSFxiHDRvmun7+8593Qf6hhx5yGBYSmKOPPto23HBDF1S\/8IUv2Fe+8hVHaCCWBGMCsF9\/qYXfddddrt8GG2xgffv2tXfeeccWLlxorB1ZCwnMpz\/9aSPIeyLVv39\/u\/XWW5uIVgjY4447zp5++mlHQrfaaiv3bPfu3Y1NobWEoBo8y+n1+9\/\/vt122202Z84cRyIJzBAW9I\/ct9xyS9nNJoRB3j\/3BKWlmOE\/Z7M56KCDnH0tWLDAOnTo4PwegkyrVcx4\/vnn7YgjjnCEdNNNNy0LP0QBv4OYPPHEEx8iMN4XIO377bef64O\/YC+f\/exno1T7k5\/8xJF2Dmaf+tSnnH9C4n7zm9+4+FTcQvEBGU444QS3eeP7v\/rVr9yB6\/rrr7eBAwdaveMD8noCgy+vWbPG\/vKXv7i\/WSPxetttt60LNr\/85S\/tnHPOMeISxMmTzfPPP9\/OPPPMqtceJXQGOjU8gSkmJMXE49xzz3Wsl02TjaU40+GddqeddrIXX3zRunTp4jbhdu3auaE8gfHj3HPPPW7zYNMuzPBAnHB6MhIEjmuuucYw4uXLlzsnghTtuOOObkxPYIqzPV72iRMn2tixY+2rX\/2qjRkzpqwZksmAnBVnYCA2ZKgIvKyFE+Qll1xSkjAQTPgcclHYPIGA4E2dOtUFH7IbbAAE3cIMB\/i99957zQgMwZkAATGDCA0YMMBhx2Y7ZMgQO+qoo8quDSwJMgRVAiEN8siJhnELCQzkivbII4\/YN77xDRekOAEdc8wxbhPi9FPq1Mozy5YtcxkiTkZsFGSXIMF\/\/OMf7eabb7YDDjggqNPChdQCz1J6vfLKKx1hYWOA6LBmCDNZGWT+3ve+l4GQlZyInqC0FDOKP8eP8WkaG5DP3NYqZrCR\/\/rXv3ZZX8g\/pJk\/HBI4MKHjAw88sBlgnNzZ+Ittmb733XefswU2zp\/+9Kd22WWXOb+bMmVK0AeYBKIDGSYDwwGATOcvfvEL57NnnHHGhxQXig\/ELtZTHIN33XVXmzlzpiNv9YwPhQSG+EUce\/nll41YDcbEBGJxKD5Ugg1ZMOKuJ0hksTlgcBDk72rXnpwX1XbmhicwsNtPfvKTTagS5El9FmdgHnvsMWe0J510kiuJ+IYRY8AQDkjAX\/\/612anlkICg0Nj8DfeeKPbjH0Z6N1333UnX5yVTAbpSU7Je+21l5EBIoBwkiFYsfl7AsPPVq1a5TJIBCCMm1YtgWGTJ+NBNuLYY491xILgBkkh+NI4RXKCopF1ImBxMqCRFeEPZTLWBCHwpR5PVAoDOs+UIjCQCdYFPjgyOLAB8yzEJJTZuPbaax354dQJoeA0s+eee7oTDa0UgeHnBHECJETpxz\/+cTBAsVERVJGHzWTrrbd2pSmIjT8xhXRaazxLERjsl\/T73\/\/+d5cVAB\/s6ZRTTnFkurCkWdswk4\/RPEFpKWa0lIHBtmfNmuUILq0WMYNx8E0IBnGFTY24BDHB1u+++27nv\/hJYQmnJQJDPCI76Qm374dNcCALbdJvvfWW8y8axJ1sE2VJSDHyEBuLWyg+QISQBxLG8z\/84Q9dzOTgw2Gx3vGhFIHhZ2TgiWfEuN\/\/\/vd1waYYqwsvvNCV\/E888UQbPXp01WvPh0eaNTyBKVYkpwJIRPEdGPoRgCA4nMppsGQ2dU6vZBRwMjZ4yglkQGilxuG0jiP6DYPNg5LKM888YwQCyAClE1KmbDSchH7+85+7MRmbuyAQKk5ZBB02YpyacgUBjEwHtVP+7YNKSwZbaqPjVM4pjJMmQYN2+OGHu4D48MMP28c+9jGHAfKWav7ExBiM5f9PX4Ip9XGPs3++FIHhM06WlDoIFpR4IBWktbmvwR8wJBBCOosb2SsIBDj6xoZNZoWA3RKB8QQRPRE00AmkiWfJIJVqBBUIlm\/IDaGBMHD6DOm01ni2lFlDPnCB\/CIjemVDpZ7PhkOj1ECJU605Ai3dgfG23NLnEA3iBifqWsUMfIKyKBldMrRs8mTU2PB79+7t4gVx48gjj3SHK99aIjD4OnekfKaBuIKPEPO4nxfyAe5sUDbDpshS0DypIjODzxa3UHygFIt9Yq8+FvhDC8SBmFDP+MCcvoTkceFnZDAhEjSwQb5y8aESbAqxIhNPDMFHIcJ+36hm7Xnx7YYnMLElJE4gvmRByYVTPZs79wrIlJBp4DIqF3rZHLgPAanwBIYNunPnzq4PdzI4XeAcNDZgggekBsdsqTRy+umnuzLN2rVrm0pU\/JvTPicSSjz+glusgZba6Mg6cMGzVIOgQdo4lUFEaGze3BPycnNRlj\/+hAUGnCBoZIvIGsVkYEoFPE6\/1OjJ9oAhd4dIW+PYLTUyIeAOGaEeTyDg\/y0RGMpTjOczMLFYYhfYCQGGC9Vki8ggcWE4pNNa41mOwPj13H777Y5wUy7Yf\/\/9XeYJe4KI\/\/a3vzUurqv9D4HWlpDIqnHg4N6EL9HWKmZwMICgYGds6v4AUKwvYgx2GCIwPgODP2ELxRmYkB1Uk4EpFx+Ii5ShIAk777yzy+a2b9\/elWp9md7LBiGqZXxoicAUZ2DqgY0fE3JMtoesLntM8Qsk1a49JHvaPxeBKbqU6xVWfHeFTZdg8N\/\/\/rcpe+CdvpSSb7rpJncqKh7H1zILAwtpYE69BDuCHqd1TvsQmsIyCRsjBICTFndoaLBwTssQGAIa9VlOBc8995zrw2msXCu10bF5s4mTveFv3yglUdPnHlBha+nOhnd0fwLzd2AgeAResku+tZSB8Z+zwYMlZTQcmY0Xske5DBLJibG4cVL6+Mc\/3lRH9qdKCCTkpRSBYXzmKbwDQ3Dm55z4iu8U+DlXr17ddOGXYMvl5H\/+858ulU\/GKqTTWuMZIjCQT+7msLFCtLkkyabIujn1Incoe5f24FZr+VpLYJgffyRbCpGl1FOrmMG9EmzOZ1ewNzIf3FvD3nxDtxzSQgTGl029nP4OjL9zEeMD3BHhbT+fkSKmkRX2ZVRiErGJww6XYkPxAaJC\/PMxkJjKm34cCJCvsNUjPpQiMKXuwNQDG+bmJQ+y7VxxILvLoblUq2TttfaNpMYTgQkQGO6gsGn7NwrY\/B5\/\/HFHNHw5oXBTg1zwh6wI9y+KCQwODmGBDJCWJE1K8GGz7dq1q7MDMi2UPXBUyllshGQOyHTwM05InDQIAqQnCQyUnAgWvOkUcweGubmgS9aAtTAPz\/KaIn+zNkpGkCoch7s9ZE3ILPlXoL3RtkRguLvCWpGfzZLgRQq98C0kLtiymVKrpzEWGFMuI4DRfJYJWdiYf\/azn7mL1Wy4nGghJqRZi9vFF1\/ssjScWlgDpzayMeDHSbOQwEAo0Q1lKhrBmCDJWkP1f\/qTbSIThA4hO9wzQa\/oIkanMQQmBs9yei0ss0F2kZdsG3pmo8K+IDDoBBtvzWuiSQWwtpw3lsD4t5DwK\/yYDZhyCr5ei5iBbUKIsW\/IJs1fDObQwp0T7rGRRaEv5WlkoXSKL3J3hAMSGUYIPpkL3mjiOXwPu\/VfB4APEHNifMDjgw9wJwyf5LDCfLwl5cusntCE7BmZIdbcW8NviW\/EBAhS4b2eesWHQgIDrhwWyab5r8Mgs0VcqQc2vJVJtonGHuTf4kJfhS9mVLr2tvSbes4lAhMgMB58HJ63gAgGZBT8PQk2ZurOvvnUK4GAOyJkMLjjUngPxG+sbBqUY6hjkyb1jRM8FywJAKRQSeVzH4PyEZdZ2RS5VMubQqRRCRbI5Z06hsD4E1excfkSEScNTjuQHOrtzEsgY+MrvjzbEoFhbEgDz+DkYMjtecpJnpz4U1uxHDznMzTUz\/kDYaEmzybLGFx0ZpNFTghYcYMU8hz3gcCRYMoaSMmSzSr1PTAEKjYZLrb6t6RiAhSBHjLAPARrNhZOs\/7V2ZBOYwhMDJ4hvTIGhI83UcDPf08N8vvX+iGH\/q2tegafrI0dS2D8urBxDgNkYfhTq5iB7ZN1Jb7416fZyMhMYu\/EBWydkguXx\/3bMxCT4laYCSYG+e+BgYT4rDDPxPiALz\/674Fh\/u985ztNmbxiAhOyZ+KOzyiTAecgSZxjPYWtXvGhkMDwbw6J4MKVAbJp\/i5KPbDx9xCL9UWmunCvqHTtWfO9luRtOALDxkwal+8B8d9FQgqQoL148WJ3T4W0HdmNtmhkFHD04leQ22JuzSEEGgUByhVk3SCSfNUBGTOyGLQk\/b8S\/LlszeGm8G5LJePoGSGQdQQaisCQ+uOVX04nlBU8geFeCqclTuYwX05avA7bFo1gRAnKf8dLW8ypOYRAoyFASp5MJadmyiMcWCiHksFM0v8r0QN31MhmUgpVEwKNjEBDERjSbZy2uH9C8wSGS268ocJryaRBeT2R9Lq\/k9LIBqK1C4G8IYD\/c5GU8ihEQP6fNw1rPY2CQMMQGL7Ei3sopJL9V+5DYLjJT7mIewH+XgZ3OrgoVe53eDSKgWidQiBvCPD6OFkMYoH8P2\/a1XoaCYGGITDc6Oa7Q3g1lItqPgPjv8m2kMBwR4a3YEKvIDeSoWitQiAPCHCxm0v0ZGO5mC7\/z4NWtYZGRaBhCAxv8fBmEM3\/dlFey+XeSyiFXPjbphvVULTubCPgfx9VtldRnfTceSHryndq+O9RYkT5f3W46un0I5BX\/28YAlNoYoUZGH7OGwm8susv8frfDO2fIdj5X\/ZXylT5MrDCL3xLvzm3TsK8rw808r7GkA23ziKy15tXXXktmK8YKP61E\/L\/+ugz7z5VKWpJ4JJX\/xeB+eA1SgiIf40aglMY5PKq\/EodUM9lD4FGt2G+T4TfV8a3KfvGHRgu7nOxX\/6fPZuWxPEI5NX\/G5LAxKv9\/3rmVfmtxUH9s4uAbLhy3Qm7yrHTk+lAIK82LAITYV95VX7E0tUlJwjIhitXpLCrHDs9mQ4E8mrDIjAR9pVX5UcsXV1ygoBsuHJFCrvKsdOT6UAgrzYsAhNhX3lVfsTS1SUnCMiGK1eksKscOz2ZDgTyasMiMBH2lVflRyxdXXKCgGy4ckUKu8qx05PpQCCvNiwCE2FfeVV+xNLVJScIyIYrV6Swqxw7PZkOBPJqwyIwEfYVUn7oi+7y+iVCEdCpS0oQCNlwSsRMpRjCLpVqkVCtQCCvNiwCE2EEIeXz+RFX39\/iSDOHDyj7RXgRIqiLEKgKgZANVzV4zh8WdjlXcAMsL682LAITYbwh5YvARICoLokiELLhRIVL+eTCLuUKknhBBPJqwyIwQdWHv8hOBCYCRHVJFIG8BrC2ALUa7MqVl1VabgvtaQ4QqMaG04ygCEyEdkLKF4GJAFFdEkUgZMOJCpfyyavBrlxsUGk55YrPkXjV2HCaYRCBidBOSPkiMBEgqkuiCIRsOFHhUj55NdiJwKRcuQ0iXjU2nGaIRGAitBNSvghMBIjqkigCIRtOVLiUT14NdiIwKVdug4hXjQ2nGSIRmAjthJQvAhMBorokikDIhhMVLuWTV4OdCEzKldsg4lVjw2mGSAQmQjsh5YvARICoLokiELLhRIVL+eTVYCcCk3LlNoh41dhwmiESgYnQTkj5IjARIKpLogiEbDhR4VI+eTXYicCkXLkNIl41NpxmiERgIrQTUr4ITASI6pIoAiEbTlS4lE9eDXYiMClXboOIV40NpxkiEZgI7YSULwITAaK6JIpAyIYTFS7lk1eDnQhMypXbIOJVY8NphkgEJkI7IeWLwESAqC6JIhCy4USFS\/nk1WAnApNy5TaIeNXYcJohEoGJ0E5I+SIwESCqS6IIhGw4UeFSPnk12InApFy5DSJeNTacZohEYCK0E1K+CEwEiOqSKAIhG05UuJRPXg12IjApV26DiFeNDacZIhGYCO2ElC8CEwGiuiSKQMiGExUu5ZNXg50ITMqV2yDiVWPDaYZIBCZCOyHli8BEgKguiSIQsuFEhUv55NVgJwKTcuU2iHjV2HCaIWoYAnPnnXfapEmTbPny5dalSxcbOXKk9e\/f3+nmrbfesqFDh9rixYutc+fONn78eOvevXuT3kLKF4FJs4lLNhAI2XDeUZo\/f76NHj3aXnjhBRs8eLANHz68acnV+n857ERg8m5Z2VhfXv2\/YQjMQw89ZLvuuqttscUW9uijjzrC8vTTT1u7du1s1KhR1r59exsxYoTNnj3bJk+ebNOnTxeByYZvSsoIBPIawCKW7rosXbrUli1bZnPnzrU1a9Y0IzDV+r8ITKwW1C8pBPLq\/w1DYAoN57333rNevXrZs88+ax07drR+\/frZtGnTrFu3brZ27Vrr27evzZkzx7p27eoeCylfGZik3FLzxiIQsuHYcbLe79prr7WVK1c2IzDV+r8ITNatIv\/y59X\/G5LA3H777TZr1iyjrLR69WpXLlq4cKHLwtAGDhxoY8aMsT59+ojA5N+3G2KFeQ1grVVeMYGphf+LwLRWC+rf1gjk1f8bjsA88cQTdvHFF9ttt91m2267ra1atcp69OjRjMAMGDDAxo4da7179xaBaWtP03x1QSCvAay1YBUTmFr4vwhMa7Wg\/m2NQF79v6EIDHdehg0bZlOmTHFlId9iUshDhgwx\/pRqKiG1tTtqvtYgMGHCBOPPkiVLWvNYLvtWWkIq5\/8iMLk0ldwsKs\/+3zAEZsGCBXb22WfbxIkTrWfPns2MkzeSOnTo0HSJl7eVZsyY0dQnxF5FYHLj67ldSMiGc7vwooWVIjDV+r8ITKNYT3bXmVf\/bxgCM2jQIJs3b5516tSpyQq5A8PFXV6j5ITlX6MeN25cM5ITUr4ITHYdu1EkD9lw3nFYtGiREQNWrFjhlrrRRhvZjTfeaLvttlvV\/i8Ck3fryf768ur\/DUNgqjHBkPJFYKpBV8+2BQIhG24LGbI6RzXYhb4HpiVM9thjD5s6dWpWIZPcKUOgGhtO2VKaiSMCE6GdkPJFYCJAVJdEEQjZcKLCpXzyarALEZgjrr6\/5OpnDh+gO0spt4ssiVeNDad5nSIwEdoJKV8EJgJEdUkUgZANJypcyievBjsRmJQrt0HEq8aG0wyRCEyEdkLKF4GJAFFdEkUgZMOJCpfyyavBTgQm5cptEPGqseE0QyQCE6GdkPJFYCJAVJdEEQjZcKLCpXzyarATgUm5chtEvGpsOM0QicBEaCekfBGYCBDVJVEEQjacqHApn7wa7ERgUq7cBhGvGhtOM0QiMBHaCSlfBCYCRHVJFIGQDScqXMonD2F3\/PHH25NPPtniKspd1NUl3pQrPyfihWw4q8sUgYnQXEj5IjARIKpLogiEbDhR4VI+eQi7arIsIjApV35OxAvZcFaXKQITobmQ8kVgIkBUl0QRCNlwosKlfPIQdiIwKVegxHO\/OiePv0pEBCbCuEPKF4GJAFFdEkUgZMOJCpfyyUPYicCkXIESTwSmkW2gmgAGbvpSqka2nnSsPWTD6ZAynVKEsBOBSafeJNX\/EAjZcFaxUgYmQnMh5SsDEwGiuiSKQMiGExUu5ZOHsBOBSbkCJZ4yMI1sA9UEMGVgGtly0rP2kA2nR9L0SRLCTgQmfTqTRM0RCNlwVvFSBiZCcyHlKwMTAaK6JIpAyIYTFS7lk4ewE4FJuQIlnjIwjWwD1QQwZWAa2XLSs\/aQDadH0vRJEsJOBCZ9OpNEysDIBj5AoJoAJgIjM0oDAiEbToOMaZUhhJ0ITFo1J7k8AiEbzipSKiFFaC6kfJWQIkBUl0QRCNlwosKlfPIQdiIwKVegxFMJqZFtoJoApgxMI1tOetYesuH0SJo+SULYicCkT2eSSCUk2YBKSLKBnCAQ2oRzssy6LCOEnQhMXWDXoDVEIGTDNZyqTYdSCSkC7pDyVUKKAFFdEkUgZMOJCpfyyUPY1YvAVALLHnvsYVOnTq3kUT2TYwRCNpzVpYvARGgupHwRmAgQ1SVRBEI2nKhwKZ88hF29CIx+i3XKDSND4oVsOENLaSaqCEyE5kLKF4GJAFFdEkUgZMOJCpfyyUPYicCkXIEST5d4G9kGqglg4MbvQgo1pX5DCOnzahAI2XA1Y+f92RB2eSAwxx9\/vD355JMlVanYlH0LD9lwVleoDEyE5kLKj8nAtJQO9tPrFz5GKEJdKkYgZMMVD9wAD4awywOBCa1hyZIlDaDp\/C4xZMNZXbkIjJm99dZbNnToUFu8eLF17tzZxo8fb927d2\/SaUj5IjBZNf\/GkTtkw42DxIdXWk\/\/52BS6V2WSp+rhGyIwOTbA\/Lq\/yIwZjZq1Chr3769jRgxwmbPnm2TJ0+26dOntzmBKedCSuPmO8DUe3V5DWC1wK2e\/p93AlNp6anS52qh76yOUQ1mefV\/ERgz69evn02bNs26detma9eutb59+9qcOXOsa9euztZDyq9VBqZcmUklpqyGnXTIHbLhdEiZjBT19P+8E5hKMzeVPpeMhaRj1mowy6v\/NzyBWb16tSsXLVy40GVhaAMHDrQxY8ZYnz59MkVgyjF0FhKTxanFGOlw99ZJMWHCBBsyZEjrHmrD3tXqJa8BrFoV1Nv\/RWBK352pZjOuVudJP19prKkGs7z6f8MTmFWrVlmPHj2aEZgBAwbY2LFjrXfv3s7WQ5tH0g6h+YVACIEY8hoaI4+fy\/\/zqFWtqRiBvPp\/wxMYFB1KIcsdhIAQyC8C8v\/86lYryzcCIjBmNnLkSOvQoUPTJd5JkybZjBkz8q15rU4ICAGHgPxfhiAEsomACMwHr1Fz\/8G\/Rj1u3Djr2bNnNjUqqYWAEGgVArxGLf9vFWTqLARSgYAITBVqIEtzzTXX2Pvvv29HHHGEXXDBBVWMlsyjb775pss8Pffcc7buuuvafvvtZ5dddlnTheaW1hj67oxkVhOe9aKLLrK77rrLFi1aZOutt557IG9rDKOgHrVAIA\/+Xw0OxL3LL7\/cvbFJ7LjwwgvtkEMO+dCQr776qu27777WpUsX9xnZ7kceeaSaqTP57PLly+2ss85ysZZLtffee28m15EmoUVgKtQGGz\/Oet999znHPOqoo5wD77PPPhWOmMxjrOP55593AWblypV26qmn2kEHHWQnn3yylVtj6LszkllN+Vkfe+wxFzT44wlM3taYRtzzKFNe\/L8a3UBcKLfz269fe+01O\/roo+2hhx6yjTfeuNmwEJiTTjrJHnzwwWqmy\/yzxNf58+fbO++8YzfccIMITA00KgJTIYg47bx581wGhsaX37EpXnXVVRWOmI7HKJ\/haGRlyq0xaxcf3333Xfc22U033WT9+\/dvIjB5WmM6LBYOzdcAACAASURBVKgxpMir\/7dGe+eee67ttddedtxxx7nHTjvtNPcVFEceeaQITBkg+Z1TV199tQhMa4ythb4iMBWCyEZPCvW8885zI3AaueOOO+yWW26pcMTkH1uxYoUdeuih7hVyXrtraY2cukLfnZP8appLMHr0aPdaPAEW2X0GJk9rTBvmeZYnj\/7fWn0de+yxds455zRlncnKbrnllnbmmWd+iMDsv\/\/+ttVWW9mGG27osjEnnHBCa6fLTX8RmNqpUgSmQiwJYLRhw4a5v++\/\/353t2LKlCkVjpjsY3wfxhlnnGG77rqr+71QtJbWOHHixOB35yS7muazk7a9\/vrrXbrbf3FZIYEppcesrTFNeDeCLHnz\/0p0dswxx7jLz3vvvbd7nLe5tt56axdHChsZ3WXLlrnPXn75ZRs0aJBdeeWVmSu3V4JRqWdEYGqFpJkITIVY5imFvGbNGiMdzAnp0ksvbUIkL+UVyAh\/1l9\/fbc2avLbbLON+\/URDz\/8cIulwKyVySo0ZT1WAQJ58v8Klu8eiS0hFY\/Pt5x36tTJZW8asYnA1E7rIjAVYukv8c2cOdM222wzd4Ft+PDhmTtV8LufyCJtsMEG7lS0zjrrNCFSbo1Z\/e6M4gxMHtdYoUnrsVYgkBf\/b8WSP9SVrDN3\/\/wlXl5k4BIv5GTBggXuMi9v2yxdutQ22WQT69ixo73++uuufMTbgAceeGA102f2WRGY2qlOBKYKLO+55x679tpr3S+APPzww91bSFlrzzzzjJEK3nzzza1du3ZO\/IMPPti4M0JraY1Z\/e6MYgKTxzVmzQazKm8e\/L8a7LkDSJyYO3euix2QEv8aNXcDd9hhB3cfhvuBV1xxhVGm5hVqLtMXl5mqkSNLz\/KGJ+U03kTadNNN7fTTT7fBgwdnaQmpklUEJlXqkDBCQAgIASEgBIRADAIiMDEoqY8QEAJCQAgIASGQKgREYFKlDgkjBISAEBACQkAIxCAgAhODkvoIASEgBISAEBACqUJABCZV6pAwQkAICAEhIASEQAwCIjAxKKmPEBACQkAICAEhkCoERGBSpQ4JIwSEgBAQAkJACMQgIAITg5L6CAEhIASEgBAQAqlCQAQmVeqQMEJACAgBISAEhEAMAiIwMSipjxAQAkJACAgBIZAqBERgUqUOCSMEhIAQEAJCQAjEICACE4OS+ggBISAEhIAQEAKpQkAEJlXqkDBCQAgIASEgBIRADAIiMDEoqY8QEAJCQAgIASGQKgREYFKlDgkjBISAEBACQkAIxCAgAhODkvoIASEgBISAEBACqUJABCZV6pAwQkAICAEhIASEQAwCDUVgTjvtNHvggQc+hMu9995ru+22m7377ru2++6723\/+8x\/beOON7amnnrIOHTo09V+4cKFdeeWV9vvf\/95WrVplW2+9tfXt29e+973vWUtj+4cfe+wxu+yyy5rmf\/\/9922LLbawfv362ciRI92\/feMzfn7ppZfaYYcdVlaP48aNsxtuuKFZn80339zJ3pp2zz332E033WR\/\/etf7SMf+Yjtvffedv3117dmiGDfF154wUaNGmXz58+3Tp062Ze+9CW74IILbL311rM\/\/\/nPdvXVV9sf\/vAH+9e\/\/mWsgc+HDh1q6667rntm+PDh9re\/\/c0+85nP2Pjx45swO\/nkk61Hjx42YsSIFmUA0x\/96Ed255132htvvOH0u8MOO9g3v\/lNt1YaNvDvf\/\/b\/XujjTayLbfc0vbZZx8766yz3L9j22233WY\/+clPbOnSpbb++uu7eVjHvvvuGztEVL+3337b4fnggw\/amjVr7HOf+5z7P7a0cuVK+853vmPz5s2z1157zeHNOi+66CL3+T\/+8Q\/71re+5T7fbrvt7Lvf\/a716dPHzQtOc+fONfxinXXWiZIlr51CMeN3v\/udffWrX3XLx1ZuueWWZlD88pe\/tAkTJtiSJUucLXz0ox+1o446yk455ZRm9laMH\/0eeeSRZn3wE+zw0EMPtfPPP9\/5hW8vvfSSHXzwwc7vC2NJKb1MnjzZpk2bZosWLXIfI19hnMFX+Nkdd9xh2FjPnj1dLPL2EdJ1W9heOVxr6cv4wqxZs+zvf\/+784Wdd97ZLrzwQqeXWrZyvsw85XSGvvB15Nxggw3s29\/+tn3lK19x4hE3iY\/se62JYbVcWz3HakgCA\/HYZpttmnBl4\/zEJz5hv\/jFL9yG5tuPf\/xjFxRoa9eudRvEq6++6jakrbbaym2mbAR\/+tOf3Abw9NNPu77PPfecvffee65f586dmzYFNlgMqVu3bi6QPfvss7ZixQo75JBD7Lrrrmua9\/nnn7cjjjjCbS6bbrppWf17AoNDsUnRNtlkE\/vBD34QbTe33367XXLJJdaxY0cbOHCgtWvXzhGKu+66K3oM35GADwk89dRTmz1LUAO\/N9980wV6Ajob\/Nlnn23Dhg0zAtI555xjn\/70px15YFOmEajPPPNMO\/roo23DDTc01vuFL3zBOSiEBjzZlOnv119KaNZCPxwc0vnOO+8YhBRn97J6AtO9e3dHCP7yl7+4v9kQZsyYYdtuu20UHhBS8CNgYBtsFGxebC7lZCw1+A9\/+EOD+BXah+93+umn269+9Svbaaed3Log1mDPWpcvX269e\/e2XXbZxdkbmyFr9pvs97\/\/fYNozZkzx5FIiBuE5fXXX3c2z0Ycu2FFgZLRTp7AtBQz2NjBkQah+O1vf+vINw0b\/+IXv2irV692umjfvr2zC+xs0qRJbmPhsESsIGbQ6Mc4Xbt2dTr3NsnfHKb8weTyyy+3E044oQlVxsNGZ8+eHUQamyeOEb\/QdzGBgXwzPnb\/qU99yvkWhP83v\/mNiy2hVm\/bC+GKfLXyZeIhvsF4HK4gMsRk\/Il42ZpWqS8zRzmdPfroo86WsEP2IA6eEBfiJYdA9hfsOI+tIQnMueee607Exa1wQ3jxxRfdZu6zEGy2BH+CE6cufzIlIEFUChsnXQIEGY0DDzyw6SMfDP38ZD3YPCBPhZkh5vz1r39t06dPd5vONddc4zZ4AgOORNDccccd3biewNC31IYzceJEGzt2rDsljhkz5kNrJrjusccetmzZMrv\/\/vtdJqO4EeiuuOIKR7ggN2QSIDylNmMw43PIRWHzBIK5pk6d6jZ1gjtBAKL2\/\/7f\/3Nje5JAJoFN1G+4EDo+g5hBhAYMGOCwY7MdMmSIO9WWa+edd54L8IWBH\/JItsfP6YMe8iHnyy+\/bIMGDXJEi0AA1gsWLHD\/3n777ZtIVrl5fTCnDzpG1yGdFo43evRoF4wgF4UN2Q444ABH9ghabI577bWXI4isk1Mz5MnbBESFTBL9fCaRMZ944gkXHNkU2PzIynDSJ6uoZk2Z1VIxo9B3IJHEDPR14oknOui8fxceUDgIkeXEfnxjU8SmaWyShRtjsU16Oz7++OOdT\/qGnUI2+DzWvtj02PyKCcxnP\/tZR2zIwPTv398d6jjc4W9nnHGGhWIKh5V62l4MrjG+DHZkVP773\/+6mF4qc8XPyYYT7zmsggcZDx9v+fett97q\/hAnOCzgT4Vx3+uoUl+G1PpWSmcPPfSQfeMb33Dxm1hKNgY7wp\/RFb6P3+exNSSB4WTMiYpGsL777rsdy+b0ivOxAXDi59SMQRBQ2OwoXfA5Gy+GjKN7IlFoHDEEhqDAhnjjjTe6zbiwDHTssce6gEYfNh0MkM2J0\/SUKVPcaQiCg1yewPAzyloQEByIjZ8WCjac7gmwPE8Q5iT48Y9\/3KVJISKsm4wHhAy5OC3ed9997jNOajQ2TXChfe1rX3PYkDWhQXL4Q1YCJ8fRfKmH+RiPE15hQOc55ifNzWaA43PiYV2U+SAc4MAGzLPoK1TquPbaax35IfNGupxMz5577ukIgG\/FQY+fU3JCXtZAhiOWwFAyfPzxxx354NT8+c9\/3umCFtIppIc\/NDaXP\/7xj44M0ziFE2gJVmSv2LTQB80HN8qcPoXs10Yf8INsk7Uh2A4ePNidKBkTfLAnShsQrVAZIo\/BsNSa\/KGjVMzAJsGQTQtfhTz4DBhjPfnkkwbRIDsG6cXm8MtibGMJTK9evZz\/YE+UoyEtNHwCoop\/MX\/Ivspthm+99ZaTk4bdUU6mjAKhxW+wx1BMKcax1rYXg2uMLyNniMAUroVSLDGfgxaZNmyCQxaHLeIX+LCXECuJSeik1r5cisBApNE5GTIatkg\/DjhkWmtduk6T7zckgSlUABsfWRTP6v2GwKbO5k7248gjj3SPsGmQ\/fjnP\/\/ZNATEgs0FR\/ctRGAK54cUQV58UINIUeJAHjZbyACpQDYdNhoC189\/\/nOXVTnuuOPcXRA2S4IoJ2tOPhAvZMWpyHTAxvm3D0yF8\/OsD4Q4PQQI4sAYbGRs2JR2OIWQAqUdfvjhzkkffvhh+9jHPmbHHHOMPfPMMyXt2p9cGYMTQeFJFgen9kt2hsDrG9kmTnoECOq6HhtOO+ADmSCYkH0hbcp9Df7QD6JE9qG4EUgoR4Gjb2zYlA8hkLRSQY8MhT9Rc8+AcSBNPAuRbakV3k1Cf5SvSPlD9kI6LXWvyc+DHiBVvuzn\/8\/nfuOCqBSWQsnWYCtkmyB+BGHfWA\/kFWzRKxstuLNp0SBJ+EKjtlJ3YHzM8NkQSAU44ctsJviUz+qhy5tvvrmJ4IMjtk02w7cYAlOIP9lUNk0OXzT8lDIsPggBCdmXH6vUZkgsPOigg5w9YDc0v0ljN\/hbKKYUylov2wvhGuPLlOogV+gMGy+M4cX2zkEOck+2s9C\/yLSQvUTHZEqIzcRo9gz2jlr7cktZM+SFyGITZFs49L3yyitOVmyFuEN2HDvJU2tIAlMqHcwmTuAhhc5m60\/s+++\/v6tX+0YKmCwFp2vuyJCuLT7xhggMmQfuxnCS4HRGqtZfCmOTxwE4ZUEeWiqNUO6iTIM8nAho\/JsLfqSyKfFwOgw1n4GhH6lIiBAbNhs3JImAyP2eUo0MDOyeNCtEhIbjcFL0cnPi54\/PwLAhkF2hQZbIGhVmYCAzZDzIkHGiLM7MeDnQEdmT\/fbbz2V7wBA9kfqG9LTUKJWBOySNEwubNf+nxZ7aQpj6z1kbWRv0wF0HTmW0kE4hSn7zQE7+7Uty1N8hez4Dg8y+vFQqA4MNff3rX3d2SmmyVGobmSBE6POnP\/2pYfP4AfbEiZvT5mabbRa77Fz1Ky77+sWRcUQP4OoJOFkv\/BayWngHjI2FTC52Dr60QpuPITBkbiDtixcvdj46c+ZMRzxp+BblDQ5C2FvIvvwaSm2GMRmYWAXX2\/bK4VpLXwZbSCoHQR93PQY+i1yMiffLWvoyc5QjMF4GCCY2QNkPssx9S54jrhdfa4jVZVr7icB8UALhzgOnjuIGo2WDhtVCcNgA\/MnHn3jZkAsvSYUIjCdQ\/p5H4SbEyYwTAezdn9Y57ROcCsskZBsgALB\/7lXQkJ+TBATGp5hxIAgXfQprqX6dhaUxUuIYuydzbF407piQveFv3yj9kK3yFxb9z1u6A+NLMYWnODIYhW97cWmRjNInP\/lJlykovGhdqBf0Ad6U0SA5bLwEblLbpEzJRBU3CBmlMX8ypg+ZFwgkWStaTN2cAM+8ZIFaIgOQOX\/ZkUvAkCxKNuiUjB0n5HI6LZS9NXVzSmJkB319Hn2SGeCkSakSGy\/VkJd0M\/om88ZGDJFGdrIypMVLZe\/SGtRqKVdLBMYTyFJzUZbgJM5GAnH1fgee6IgN0b\/5yPMxBIZ7WRAm9MGho3AjpdyMnsmyxcSMcgSGz7AT3tTzxIyxWa+\/UB+KKYxRT9uLwTXGl5GTww4HDeJm4RunHiN0QzzkEMHdEg4DhY3SMJ9xH4lstG\/4EdnRWvlySGeF80BUyG6TfSGWUjWAwOy6667Obgqzs7X0lSTGEoExc6ciNvwuXbo0uwjLpVVO7BgnaVU2HjYu7r1AMjhhQBoIRoXkIJbAECQIPpABf3GUOS6++GIXqGj+YjFjEsDYoMgcsLHxM07LZCJ4c4b0L2SFkgUBBzISU68mg0QKFIfjng14kCGilANxg2hQMiJAQiq4hMgJknSqJ3MhAsNpibUiP5slQZAUp38LyV9EYxxk8AQA0lF4+dhnmZAFp\/zZz37mLmSz4VLigpggd3EDU7I0ZHRYA5sAui3MsPmgB9lDr7yF5F93R8c8F3MHhnGwB0gm8\/AHokT5CgIV0mlM0KOP31zBi\/EpIUA0\/J0u\/o2d8pp0YXCFzBTe\/cH2yVyRUUTPnPaxLzZcLj2TbczjK5gxAbclAuM3dbKLhVlCdMxBA\/\/gbUJIP9lE8Mff2eyIM\/iwL1nEEhiIhb+MDQFGL5AhNlEIOm8uxcQMYg2xjQu8EB7iCrGCzZnyK6\/sEvMYj00Pf+KgASkhAxiKKZR562l7kI4QrjG+DFahOzCUgsi8sH6w8I0MJWv0b2wRM7gDw4GQe2\/oqrhc09JhJOTLfB7SmZeLuI\/uKCtiX9gGB00ygpA0srB5KgmLwPz\/l9X8HY7iTIq\/vMbGTh2YEypOzwYOaydwEeD8HRlvRLEEhv5+Y2XTYH6Mi43Ivz7NHQXKWQQRghWpfGrtlBTYDAkm1HF5U4hSEk6G4\/i3T0LBBhlIh0NiIANseDg\/cvm7JATdq666ygU9ghPzskY2vuLLsy1lYJiHjZxnIAEEYN4sYs2QJAgIZaHiRoqWVKhv1OD5Q38yC2yyjEFgZ5NFzsL7NP45AjzPEYzAkYDEGihX+Ts2hd8DAwkkgJMx4cTi+8QQGNZBKYFsDUSMIMkYXq6QTmMJDOODJxkh9IYNkdVDbsgZdlKqka3yb5BB+LiwDn4QHRoE1r+2DTksfF03ZtPPU59SBIayEdhy8KCUU7ix+RQ\/WVbS+JBCMoYQQ07lkE2yrIUHntYQGAg1hykyr2REyBpwX67QR0L2xeaL3MXNl4R96dB\/Dwz2S\/bBZ+FCMaXetkesC+Ea48sxBIZ7dpTtipv\/ig2wIgtMbEEu4gp4kZ0vzniWIzDlfJm5QzrzcRzCwt7gv9cH2yPGQZ65XsBh0F85yIOfNhSBSbvCIEwQleJXZtMut+QTAmlHAMLLhUpemSVbxQbz5S9\/2YnN5gFRY6Mi80gZk4xmFhqEiQ2z+GsLsiC7ZBQC1SIgAlMtgjV8ntQomQnKGmpCQAjUDgGfdSNjwX0kMoWUebjDRdaKLCAnVYgOJRTuEWWhcejhLlapr3PIgvySUQhUg4AITDXo6VkhIAQyhwB3myhf+u+7oUTM21683UNJgK8xoKTm75RkboESWAg0CAIiMA2iaC1TCDQ6AmRXyLbwfTjcBeDuG3eHKBdxF8h\/WynZGT7Xr1JodIvR+tOOgAhM2jUk+YSAEKgpAlwm5y0WXhOmZMvXERQSGF6v5zuQSn3tQE0F0WBCQAhUhYAITAR8fB24\/7KziO7qIgRSh4D\/HVSpEywhgbj8yttVvGUSKiHJ\/xNSkqatGQJ59X8RmAgT4R1\/XltsqfF6deGXvEUM2aZdkP+Iq+8Pzrnwgdutx4H\/+w23LT0wc\/iAsngEJ2qDDmnXSSwEtVpHyIZj5clqP74fiUwLrzLzRWh8jxAXdcGFb7Llcq+\/xOt\/s7Nfawi7cv6VBV9pC53Wyo7bQtYk56gXTiEbTnLN1cwtAhOBXtaVH0tgIqBwXRSUY5FKT7+s23C1SPJdSnwbNI3vwOE1av\/9TbxGzQHEv0bN69aFv08rhJ0ITLXa0fP1RiBkw\/Wev17ji8BEIJt15YvARCg5512ybsNJqieEnQhMktrR3DEIhGw4Zow09hGBidBK1pVfDwITAVt0l7zWZ6MBaIOOWbfhNoCoxSlC2InAJKkdzR2DQMiGY8ZIYx8RmAitZF359SAwMXdqIqB1XVSSikWq8n5Zt+HKV179kyHsRGCqx1gj1BeBkA3Xd\/b6jS4CE4Ft1pUvAhOh5Jx3yboNJ6meEHYiMElqR3PHIBCy4Zgx0thHBCZCK1lXvghMhJJz3iXrNpykekLYicAkqR3NHYNAyIZjxkhjHxGYCK1kXfkiMBFKznmXrNtwkuoJYScCk6R2NHcMAiEbjhkjjX1EYCK0knXli8BEKDnnXbJuw0mqJ4SdCEyS2tHcMQiEbDhmjDT2EYGJ0ErWlS8CE6HknHfJug0nqZ4QdiIwSWpHc8cgELLhmDHS2EcEJkIrWVe+CEyEknPeJes2nKR6QtiJwCSpHc0dg0DIhmPGSGMfEZgIrWRd+SIwEUrOeZes23CS6glhJwKTpHY0dwwCIRuOGSONfRqGwMyfP99Gjx5t\/CbawYMH2\/Dhw5v0wVeJDx06tOmrxMePH2\/du3dv+jzryheBSaPrta1MWbfhtkWr+Wwh7ERgktSO5o5BIGTDMWOksU\/DEJilS5fasmXLbO7cubZmzZpmBGbUqFHWvn37pl\/mNnnyZPeL3nzLuvJFYNLoem0rU9ZtuG3REoFJEm\/NXXsE8ur\/DUNgvEnwS91WrlzZjMD069fPpk2bZt26dbO1a9da3759bc6cOda1a1f3WNaVLwJT+4CQtRGzbsNJ4h3CThmYJLWjuWMQCNlwzBhp7NPwBGb16tWuXLRw4UKXhaENHDjQxowZY3369BGBKWG1fPW\/fpVAGt25ZZnyGsBitcBvop40aZItX77cunTpYiNHjrT+\/fu7x6stIYvAxGpB\/ZJCIK\/+3\/AEZtWqVdajR49mBGbAgAE2duxY6927twiMCExSMaem8+Y1gMWC9NBDD9muu+5qW2yxhT366KPuztvTTz9t7dq1s2pLyCIwsVpQv6QQyKv\/NzyBwaBiSkhDhgwx\/mSxqYSURa3VTuYJEyYYf5YsWVK7QTM80nvvvWe9evWyZ5991jp27Bjl\/+WwE4HJsDE0iOgiMDlRdKk7MKSTO3To0HSJl1TzjBkzmlacdeWLwOTEeKtYRtZtuIqlf+jR22+\/3WbNmmWUlWpRQhaBqaV2NFY9EMir\/zdMBmbRokU2aNAgW7FihbOPjTbayG688UbbbbfdXA2c7MrixYutc+fONm7cOOvZs6cITAuepDsw9Qgx9R0zrwGstag98cQTdvHFF9ttt91m2267rdWihCwC01otqH9bI5BX\/28YAlONwWRd+crAVKP9fDybdRuuhRa48zJs2DCbMmWKe7PQt2pLyCIwtdCOxqgXAnkuIYvARFhN1oO\/CEyEknPeJes2XK16FixYYGeffbZNnDixWXaVcastIYvAVKsdPV9vBPLq\/yIwEZaTdeWLwEQoOeddsm7D1aqH8vG8efOsU6dOTUNxB4bvfqq2hCwCU6129Hy9Ecir\/4vARFhO1pUvAhOh5Jx3yboNJ6meEHYiMElqR3PHIBCy4Zgx0thHBCZCK1lXvghMhJJz3iXrNpykekLYicAkqR3NHYNAyIZjxkhjHxGYCK1kXfkiMBFKznmXrNtwkuoJYScCk6R2NHcMAiEbjhkjjX1EYCK0knXli8BEKDnnXbJuw0mqJ4SdCEyS2tHcMQiEbDhmjDT2EYGJ0ErWlS8CE6HknHfJug0nqZ4QdiIwSWpHc8cgELLhmDHS2EcEJkIrWVe+CEyEknPeJes2nKR6QtiJwCSpHc0dg0DIhmPGSGMfEZgIrWRd+SIwEUrOeZes23CS6glhJwKTpHY0dwwCIRuOGSONfURgIrSSdeWLwEQoOeddsm7DSaonhJ0ITJLa0dwxCIRsOGaMNPYRgYnQStaVLwIToeScd8m6DSepnhB2IjBJakdzxyAQsuGYMdLYRwQmQitZV74ITISSc94l6zacpHpC2InAJKkdzR2DQMiGY8ZIYx8RmAitZF35IjARSs55l6zbcJLqCWEnApOkdjR3DAIhG44ZI419RGDM7MILL7Q5c+bYBhts4HR01llnGb87xbesK18EJo2u17YyZd2G2xat5rOFsBOBSVI7mjsGgZANx4yRxj4iMB8QmL322ssOO+ywkjrKuvKzQGBq7Rx77LGHTZ06tdbDZna8rNtwksCHsBOBSVI7mjsGgZANx4yRxj4iMCIwrbbLmcMH2BFX39\/q51p6oNbjMQ9jLlmypGYyZn2gvAawWL3Mnz\/fRo8ebS+88IINHjzYhg8f3vQov4166NChtnjxYuvcubONHz\/eunfvHp2BFYGJ1YL6JYVAXv1fBOYDAvPwww\/beuutZ7169bJLLrnEtt122+gAlpRRxs6bhQxMLQmRCMyHLSOvASzWB5YuXWrLli2zuXPn2po1a5oRmFGjRln79u1txIgRNnv2bJs8ebJNnz492v9FYGK1oH5JIZBX\/xeBMTOC25Zbbuls64YbbrBHH320VQEsKaOMnVcEJhap\/PbLawBrrcauvfZaW7lyZTMC069fP5s2bZp169bN1q5da3379nV34rp27eqGD2EnAtNaLah\/WyMQsuG2lqdW84nAFCH5zjvvuAC2cOFCW2eddZoC2JAhQ4w\/WWwiMFnUWu1knjBhgvFHJTWzYgKzevVqVy7C38nC0AYOHGhjxoyxPn36iMDUzgw1UoIIiMAkCH69p37ppZdsxx13tPfff9+lj++77z6bNWtWdAq53vJVO74ITLUIZv\/5vAaw1mqmmMCsWrXKevTo0YzADBgwwMaOHWu9e\/cWgWktwOqfSgTy6v\/KwJjZiSee6AJYu3btXDC77LLLXNrYt6wrXwQmlTGlTYXKug3XCqxKS0jlMrAqIdVKOxqnHgjkOQMrAhNhMVkP\/iIwEUrOeZes23Ct1FOKwIwcOdI6dOjQdIl30qRJNmPGjOgDjAhMrbSjceqFQF79XwQmwmKyrnwRmAgl57xL1m24WvUsWrTIfTnlihUr3FAbbbSR3XjjjbbbbrsZr1GTYfGvUY8bN8569uwpAlMt6Ho+NQjk1f9FYCJMLOvKF4GJUHKgy\/HHH29PPvlk9QN9MMLmO3zKnvrVfTUbLzRQ1m04tL56fh7CThmYTojqrQAACi5JREFUeqKvsWuBQMiGazFHEmOIwESgnnXli8BEKDnQpR4YtuVbQVm34eo1WPkIIexEYCrHVk+2DQIhG24bKWo\/iwhMBKZZV349Nt9afvFcFr6Jtx4YisBEOF8KuoT8XwQmBUqSCGURCNlwVuETgYnQXNaVX4\/NVwQmwnDKdGnrX3WQdRuuDu3qng5hJwJTHb56uv4IhGy4\/hLUZwYRmAhcs658EZgIJauEVD1IOR0h5P8iMDlVfI6WFbLhrC5VBCZCc1lXvghMhJJFYKoHKacjhPxfBCanis\/RskI2nNWlisBEaC7ryheBiVCyCEz1IOV0hJD\/i8DkVPE5WlbIhrO6VBGYCM1lXfkiMBFKFoGpHqScjhDy\/xCBaQmWPfbYw6ZOnZpT1LSsNCEQsuE0ydoaWURgItDKuvJFYCKULAJTPUg5HSHk\/yEC09KF97a+yJ1T9WhZEQiEbDhiiFR2EYGJUEtbK7\/WX5rGEtP81lC9XqOOUG2rutQaQ71G3Sr4E+sc8n8RmMRUo4kjEQjZcOQwqesmAhOhkrZWfqNlTOpFYGpNOGo9nghMhPOloEvI\/0VgUqAkiVAWgZANZxU+EZgIzbW18kVgIpQS6FJrUlSP8URgqtdzW4wQ8v96EZhymVjdn2kLzednjpANZ3WlIjARmmtr5YvARChFBKZVILW1DbdKuJR3DmFXLwITGrctCXDKVSTxAgiEbDirAIrAmLnfRjt06NCm30Y7fvx46969e5NOQ8rXnZXqzL\/W2Q2kqfWY9RivOtSaPx06kYdsuJayZG2sav0\/RDTKXeINYaULwCGE9HkMAnn1fxEYMxs1apS1b9\/eRowYYbNnz7bJkyfb9OnTowmMMiYxLtRyn1qTg6wQmLa8U5PXAFad5f3f0\/X0\/3K2HbL70LPKwNRC+40xRl79XwTGzPr162fTpk2zbt262dq1a61v3742Z84c69q1q7PukPJFYKoLAqFAXsnotR4zC+OV29BCNlwJxnl5pp7+HyIh5Uhs6NmW8A9l4\/KityTXkcT9pGrmzKv\/NzyBWb16tSsXLVy40GVhaAMHDrQxY8ZYnz59GorALHzgdutx4AnBuJD2zTxPGZjW6EQEJmi6H+pQb\/8PkZBqCEyWyksTJkywIUOGtF5BKX0iVDasNDtWDqdq5hSBSakhVSvWqlWrrEePHs0IzIABA2zs2LHWu3dvN3w97rhUK7eeFwKtQUCn8tJoyf9bY0Xqm1UE8ur\/DZ+BwSBDKeSsGq3kFgJCIIyA\/D+MkXoIgTQiIAJjZiNHjrQOHTo0XeKdNGmSzZgxI436kkxCQAjUGAH5f40B1XBCoI0QEIH54DVq6rOLFy+2zp0727hx46xnz55tpAJNIwSEQJII8Bq1\/D9JDWhuIVAZAiIwleHmniJLc80119j7779vRxxxhF1wwQVVjFb\/R\/fdd1979913rV27dm6ym2++2Xr16lV2LaHvyKi\/1Gbz58+30aNH2wsvvGCDBw+24cOHN01bTr7nn3\/ezj\/\/fHv77bdtl112Mb7fp2PHjonqrtxaLrzwQvf22wYbbOBkPOuss2zQoEHu32lcS1voPs1zZM3\/2wrLSu24reRLap7vfOc79tBDD9nrr7\/uYhlZf99asqU0xN+k8IqZVwQmBqUSfd5880075JBD7L777rMuXbrYUUcdZTjuPvvsU+GI9X8MAnPPPffYFlts0WyycmsJfUdG\/aU2W7p0qS1btszmzp1ra9asaUZgysnH22TnnXeeHXDAAXbppZdap06d3P+T1F25tWA\/e+21lx122GEfgjWNa2kL3ad1jiRtKK2YeLkqseO0r6kW8j311FO2\/fbbW\/\/+\/ZsRmLTH31qsvV5jiMBUiOzUqVNt3rx5LgND48vvFi1aZFdddVWFI9b\/sZYITLm1pOmC47XXXmsrV65sRmBako9M04knnmiPP\/64A5YTz9lnn+1OQGnQXam1tBT4X3755VSvpf6Wm74Z0mBD6UPl\/ySqxI7TupZ6yMUrzYUZmKzE33pgUe2YIjAVIsg9GUpHnOhppP7vuOMOu+WWWyocsf6PQWBolJD2339\/RwRIY7a0Fi4zh74jp\/5S\/2+G4k2\/3Hd4rFixwpWM\/DcqU0bi5PPiiy+2uN621F1LBObhhx+29dZbz5X2LrnkEtt2223tscceS\/Va2tIG0jJXFv2\/rbCDwLTWjttKtjTMU0xgshJ\/04BdsQwiMBVqBaOjDRs2zP19\/\/3321133WVTpkypcMT6P\/bKK6\/YdtttZ8uXL3dy77TTTo6AtbSWiRMnBr8jp\/5St0xgyn2HxzvvvGN8KdTdd9\/tBqCWvOeeezYRmKR1V4rAUF7acsstnbw33HCDPfroo46AQWDSvJa2tIG0zJVF\/28r7Cqx47aSLQ3zlCIwpeJR2uJvGrATgamRFrKeQn7ggQdctujWW28tW1LJcgnphBNOsCeeeMJpPAslpELThIDxKy34hui\/\/OUvlua11MilMjVM1v2\/rcCOteO2kicN86iEVDstKANTIZb+4tXMmTNts802s6OPPtqVZNJ6iZesy3\/+8x\/3+524H8LbOdtss4377ptya0nTd2SUylqUk49vVGad\/hIvbyDxplgadFdqLS+99JLtuOOOrjTJnSouiM+aNctZaJrXUqELZfqxNNhQWgGsxI7TupZ6yFVMYLISf+uBRbVjisBUgSBv9LAR8QsgDz\/8cHd5La2N8tEpp5ziykfrrruue9uFzd+\/VtzSWtLwHRlcjuZ1Yu610DbaaCO78cYbbbfddnOloZa+w2PBggWOVLLmnXfe2V245k0kWlK6K7cWLh2TceGOEr\/e4rLLLnO\/SJSWxrWk1dbbSq6kbKit1lfpPJXacaXzZeU5fzfojTfecKXi3Xff3a677rqy8SgN8TfN+IrApFk7kk0ICAEhIASEgBAoiYAIjAxDCAgBISAEhIAQyBwCIjCZU5kEFgJCQAgIASEgBERgZANCQAgIASEgBIRA5hAQgcmcyiSwEBACQkAICAEhIAIjGxACQkAICAEhIAQyh4AITOZUJoGFgBAQAkJACAgBERjZgBAQAkJACAgBIZA5BERgMqcyCSwEhIAQEAJCQAiIwMgGhIAQEAJCQAgIgcwhIAKTOZVJYCEgBISAEBACQkAERjYgBISAEBACQkAIZA4BEZjMqUwCCwEhIASEgBAQAiIwsgEhIASEgBAQAkIgcwiIwGROZRJYCAgBISAEhIAQEIGRDQgBISAEhIAQEAKZQ0AEJnMqk8BCQAgIASEgBISACIxsQAgIASEgBISAEMgcAiIwmVOZBBYCQkAICAEhIAREYGQDQkAICAEhIASEQOYQEIHJnMoksBAQAkJACAgBISACIxsQAkJACAgBISAEMoeACEzmVCaBhYAQEAJCQAgIAREY2YAQEAJCQAgIASGQOQREYDKnMgksBISAEBACQkAIiMDIBoSAEBACQkAICIHMIfD\/ATkiNwkPd89CAAAAAElFTkSuQmCC","height":337,"width":560}}
%---
