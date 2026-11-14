function results = analyzeEEG(o,ar,pv)
% Example analysis comparing EEG/LFP signals with TCS artifacts with
% signals that are recovered with artifact removal.
arguments
    o (1,1)                         artSim % The artSim object that simulated the signals and artifacts
    ar    (1,1) struct                  % The recovered signals (results output from artifactRecovery output).
    pv.freqRes (1,1) double         = 1
    pv.highCutOff (1,1) double      = 35 % Hz
    pv.lowCutOff (1,1) double       = 1 % Hz
    pv.tag (1,1) string             = ""
    pv.exportFig (1,1) logical      = false
    pv.errorThreshold (1,1) double  = 1 % More than 1% error is "forbidden".
    pv.amplitudeThreshold (1,1) double =0.1e-6 % Ignore frequencies with amplitudes less than this [V].
    pv.fillHz (1,1) double          = 1;
    pv.axs  = NaN
    pv.tlim = [0 0.1] + mean(o.tRecord);
    pv.lineWidth = 2;
    pv.legend (1,1) logical = true;
    pv.showErrorMuv (1,1) logical = false
end

%% Analyze the power spectrum
[pwr,frequency] = pspectrum([o.vTruth ar.signal o.vContaminated ],o.recordingSamplingRate,'FrequencyLimits',[pv.lowCutOff pv.highCutOff],'FrequencyResolution',pv.freqRes);
amplitude = sqrt(2)*sqrt(pwr); % Convert to V amplitude.  (pspectrum is single sided in rms,hence *sqrt(2))

% Determine error in frequency domain
results.mismatch = amplitude(:,2)-amplitude(:,1);
% Divide by neural truth to get pct error
pctError = 100*(results.mismatch)./amplitude(:,1);
% Ignore frequencies where the signal is below the amplitudeThreshold
% to avoid large pctError for empty frequency bins
hasSignal = amplitude(:,2)> pv.amplitudeThreshold;
pctError(~hasSignal) = NaN;
meanPctError = mean(abs(pctError),"omitmissing");
results.amplitude = amplitude;
results.pctError = pctError;
results.frequency = frequency;
results.meanPctError  = meanPctError;

%% Determine the forbidden zones in the spectrum
nrToFill = pv.fillHz/(results.frequency(3)-results.frequency(2)) ;
[~,~,from, to] = fillRegions(abs(results.pctError)>pv.errorThreshold,nrToFill);
forbidden = (results.frequency(to)-results.frequency(from));
results.forbidden = sum(forbidden);
if o.tacsFrequency==0
    inTacsBand = false;  % in-band undefined for tDCS.
else
    if ~isempty(from) & ~isempty(to)
        inTacsBand = from<find(results.frequency<o.tacsFrequency,1,"last") & to > find(results.frequency > o.tacsFrequency,1,"first");
    else
        inTacsBand =false;
    end
end
if any(inTacsBand)
    results.forbiddenBand = forbidden(inTacsBand);
else
    results.forbiddenBand = 0;
end


%% Generate a figure
if ~isempty(pv.axs)
    if any(isgraphics(pv.axs))
        axs = pv.axs;
    else
        [f,axs] = fig('Name',strjoin(arMode,'/')+pv.tag,'paperCols',1,'height',9,'nrRows',3,'byColumn',false);
    end

    if numel(axs)>=1 && isgraphics(axs(1))
        axes(axs(1)) % Time course
        yyaxis left  % Truth and recovered signal (small)
        tStay = o.tRecord >= pv.tlim(1)  &  o.tRecord < pv.tlim(2) ;
        y = [o.vTruth(tStay) ar.signal(tStay)]./1e-6;
        y = y - mean(y);
        x = o.tRecord(tStay)-pv.tlim(1);
        h = plot(x,y,'LineWidth',1,'LineStyle','-');
        ylabel 'Amplitude (\muV)'
        ax =gca;
        ax.YColor ='k';
        hold on
        plot(xlim,zeros(1,2),'k-')
        ylim([-1 1]*max(abs(ylim)))
        yyaxis right % Signal with artifact
        y= o.vContaminated(tStay)/1e-6;
        y= y-mean(y);
        h(3) = plot(x,y,'LineWidth',1,'LineStyle','-');
        h(1).Color = 'r';
        h(2).Color = 'g';
        h(3).Color = 'b';
        h(1).LineWidth =pv.lineWidth;
        ylim([-1 1]*max(abs(ylim)))
        if pv.legend
            legend(h,'Neural','Recovered','Recorded','Location','NorthEast')
        end
        ax =gca;
        ax.YColor ='b';
        xlabel 'Time (s)'
        ylabel 'Amplitude (\muV)'

    end
    if  numel(axs)>=2 && isgraphics(axs(2))
        axes(axs(2)) % Spectrum
        yyaxis left  % Truth and recovered signal (small
        x = frequency;
        y= amplitude(:,[1 2])./1e-6;
        h = plot(x(1:10:end),y(1:10:end,:),'LineWidth',1,'LineStyle','-');
        ylabel 'Amplitude (\muV)'
        hold on
        ax =gca;
        ax.YScale = 'Log';
        ax.YColor ='k';

        yyaxis right % Signal with artifact
        h(3) = plot(frequency,amplitude(:,3)/1e-6,'LineWidth',1,'LineStyle','-');
        h(1).Color = 'r';
        h(2).Color = 'g';
        h(3).Color = 'b';
        h(1).LineWidth =pv.lineWidth;
        plot(xlim,zeros(1,2),'k-')

        if ~isgraphics(axs(1))  && pv.legend
            legend(h,'Neural','Recovered','Recorded')
        end
        ax =gca;
        ax.YScale = 'Log';
        ax.YColor ='b';
        set(gca,'XScale','Linear','XLim',[pv.lowCutOff pv.highCutOff]);
        xlabel 'Frequency (Hz)'
        ylabel 'Amplitude (\muV)'
    end
    if  numel(axs)>=3 && isgraphics(axs(3))
        axes(axs(3)); % Error as a percentage
        if pv.showErrorMuv
            yyaxis left
        end
        plot(frequency,pctError,'LineWidth',1,'Color','k');
        xlabel 'Frequency (Hz)'
        ylabel 'Error (%)'
        set(gca,'XScale','Linear','YScale','Linear','XLim',[pv.lowCutOff pv.highCutOff],'YColor','k');
        ylim([-1 1]*max(abs(ylim)))
        hold on
        plot(xlim,zeros(1,2),'k-')
        if pv.showErrorMuv
            yyaxis right
            plot(frequency,results.mismatch/1e-6,'LineWidth',1,'Color','m')
            set(gca,'YScale','Linear','YColor','m');
            ylabel 'Error (\muV)'
            ylim([-1 1]*max(abs(ylim)))
        end
    end
    % Export figure for paper
    if pv.exportFig
        name = sprintf('%s.png',pv.tag);
        exportgraphics(f,fullfile('../docs/',name),"Colorspace",'rgb','ContentType','vector','Resolution',600);
    end

end
%  fprintf('%s (r=%.2f rmse = %.0f/%.0f muV ae = %.1f%%)\n',strjoin(arMode,'/'),results.r.(arMode),results.mse.(arMode),o.additive/1e-6,meanPctError);
