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
    for thisRestData = fetch(allRestData,'LIMIT 2')'
        sCntr = sCntr+1;
        % Extract a Fieldtrip dataset
        ft = ephys.fieldtrip.dataset(thisRestData,dimension="rest");
        % Add the data to the ArtSim object (focus on electrodes near Oz)
        readData(A,ft,format="FIELDTRIP",channel= [125 126 138 139 117]);         %[output:4d3071ed]
        for arMode = arModes
               
            r  = simPlotEEG(A,ar,arMode,freqRes =1,axs=[],legend=false); %[output:9a247f74]
            s10 = r.amplitude(r.frequency>9 & r.frequency <11,1:2);
            s5 = r.amplitude(r.frequency>4 & r.frequency <6,1:2);
            
            thisPct5 = mean(100*abs((s5(:,2)-s5(:,1))./s5(:,1)),"omitmissing");
            thisPct10 = mean(100*abs((s10(:,2)-s10(:,1))./s10(:,1)),"omitmissing");        
            pctError = [pctError; table(string(thisRestData.subject),arMode,thisPct5,thisPct10,'VariableNames',{'subject','arMode','e5','e10'})];            %#ok<AGROW> %[output:9c386a7d]
        end
    end
end %[output:group:70b5a62b]

[pctError(:,:,1),ix] = sort(pctError(:,:,1));
pctError(:,:,2)  = pctError(ix,2);
figure;
yyaxis left
ylabel '%Error'
plot(pctError(:,1))
hold on
yyaxis right 
plot(pctError(:,2))
xlabel 'Subject'
ylabel '%Error'
legend('5Hz', '10Hz')
title 'Resting State with 10Hz tACS'
%[text] ## 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:4d3071ed]
%   data: {"dataType":"text","outputData":{"text":"the call to \"ft_selectdata\" took 0 seconds\nUsing 240.0 s of data from the FIELDTRIP struct.\n","truncated":false}}
%---
%[output:9a247f74]
%   data: {"dataType":"text","outputData":{"text":"RBAR (r=1.00 rmse = 0\/0 muV ae = 3.8%)\n","truncated":false}}
%---
%[output:9c386a7d]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('table', 'C:\\Program Files\\MATLAB\\R2025a\\toolbox\\matlab\\datatypes\\tabular\\@table\\table.m', 87)\" style=\"font-weight:bold\">table<\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025a\\toolbox\\matlab\\datatypes\\tabular\\@table\\table.m',87,0)\">line 87<\/a>)\nIncorrect number of arguments.\n\nCaused by:\n    You might have intended to create a one-row table with the character vector '288' as one of its variables. To store text data in a table, use a string array or a cell array of character vectors rather than character arrays. Alternatively, create a cell array with one row, and convert that to a table using CELL2TABLE."}}
%---
