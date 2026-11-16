function [sortResults,falsePositives,falseNegatives,spikes] = analyzeSpikes(o,ar,pv)
            % Helper function to run the simulation, recording and artifact correction
            % and show the results for an single unit recording experiment. See
            % spikeArtifacts.mlx for examples.
            arguments
                o (1,1) artSim
                ar (1,1) struct       % Artifact removal parms                
                pv.lowCut (1,1) double = 300
                pv.highCut (1,1) double = 3000 % Hz
                pv.filterOrder (1,1) double = 1
                pv.threshold (1,:) double =1
                pv.sharedThreshold (1,1) logical= true
                pv.spikeSortingMode (1,1) string = "MATCHSHAM"
                pv.nrPhaseBins (1,1) double = 8
                pv.phaseFromRaw (1,1) logical = true;
                pv.showMissed (1,1) logical  = false
                pv.tag (1,1) string = ""
                pv.exportFig (1,1) logical = false
                pv.axs  = []
                pv.xlim  = 4+[0 0.250];
                pv.nrBootstraps (1,1) double  =0
            end
            colors = 'br';
          
            

            %% Sort spikes
            [sortResults,spikes] =sortSpikes(o,ar.signal,...
                'sortMode',pv.spikeSortingMode, ...
                'spikeBand',[pv.lowCut pv.highCut], ...
                'filterOrder',pv.filterOrder, ...
                'threshold',pv.threshold, ...
                'sharedThreshold',pv.sharedThreshold, ...
                'phaseFromRaw',pv.phaseFromRaw, ...
                'nrPhaseBins',pv.nrPhaseBins,...
                'nrBootstraps',pv.nrBootstraps);

            %% Visualize
            % Generate a figure
            arMode= strjoin(cellfun(@(x)(x.name),ar.step),'/');
            if isempty(pv.axs)
                [f,axs] = fig('Name',"Spikes" + arMode,'paperCols',2,'height',9,'nrRows',2,'nrCols',2,'byColumn',true);
            else
                axs = pv.axs;
            end

            % A/B - Filtered voltage and spike detection
            for i=1:2
                if isgraphics(axs(i))
                    axes(axs(i)); %#ok<LAXES>
                    h0= plot(o.tRecord,sortResults(i).voltage/(o.amplifierGain*1e-6),'LineWidth',0.5);
                    h1 = plot(o.tSimulate(o.ixSpike),zeros(size(o.ixSpike)),'go','MarkerSize',7); % Ground truth
                    h2= plot(o.tRecord(sortResults(i).ixDetected),zeros(size(sortResults(i).ixDetected)),'g*','MarkerSize',7);
                    ylim(o.spikePeak/1e-6*[-1 1])
                    plot(pv.xlim,spikes.info.detect.thresh(min(i,numel(spikes.info.detect.thresh)))/(o.amplifierGain*1e-6)*[1 1],'k')
                    xlabel 'Time (s)'
                    ylabel ('Filtered Voltage (\muV)')
                    if i==1
                        legend([h0 h1 h2],'V','Spk','Detect','Orientation','horizontal','Location','north')
                    else
                        h3= plot(o.tSimulate,max(ylim)*o.vTcsActual/max(o.vTcsActual),'b:','LineWidth',1);
                        legend (h3,'tACS')
                    end
                    xlim(axs(i),pv.xlim)
                end
            end
            % C  - Polar spike histogram
            if numel(axs)>2 && isgraphics(axs(3))
                axs(3).Visible = 'off';
                axs(3) = polaraxes('Position',axs(3).Position);

                h = nan(1,2);
                maxRho =0;
                for i=1:2
                    theta   = [sortResults(i).binCenters;sortResults(i).binCenters(1)];
                    rho     = [sortResults(i).spikeCount;sortResults(i).spikeCount(1)]./[sortResults(i).phaseCount;sortResults(i).phaseCount(1)];
                    h(i)= polarplot(axs(3),theta,rho,'Color',colors(i),'LineStyle','-','LineWidth',1);
                    hold on
                    %  polarplot(axs(3),theta,rho, '.','Color',colors(i));
                    maxRho = max([maxRho;rho]);
                end
                axs(3).ThetaTick = [0 45 75 90 135 180 225 270 315];
                axs(3).ThetaTickLabel = {'0','45','','90','135','180','225','270','315'};
                axs(3).RLim = [0 maxRho];
                axs(3).RTick = [0 maxRho];
                axs(3).RTickLabel = [];
                axs(3).RAxisLocation =75;
                legend(h, 'Sham','tACS','Location','SouthEastOutside')
            end

            % D - waveforms for each phase
            if numel(axs)>3 && isgraphics(axs(4))
                axes(axs(4))
                axs(4).Visible = 'off';
                if pv.showMissed
                    m =   cat(3,sortResults(1).meanMissedWF(25:45,:),sortResults(2).meanMissedWF(25:45,:));
                    sd = cat(3,sortResults(1).sdMissedWF(25:45,:),sortResults(2).sdMissedWF(25:45,:));
                else
                    m =   cat(3,sortResults(1).meanWF(25:45,:),sortResults(2).meanWF(25:45,:));
                    sd = cat(3,sortResults(1).sdWF(25:45,:),sortResults(2).sdWF(25:45,:));
                end

                plotWaveformPhase(sortResults(1).binCenters, m,sd,...
                    'ax',axs(4),'UseBins',1:2:pv.nrPhaseBins, ...
                    'radius',0.25,'width',0.4,'height',0.4, ...
                    'colors',colors);
            end

            if pv.exportFig
                name = sprintf('%s.png',pv.tag);
                exportgraphics(f,fullfile('../docs/',name),"Colorspace",'rgb','ContentType','vector','Resolution',600);
            end

            % Summarize results on the command line
            trueNrSpikes = numel(sortResults(1).found); % Same for both
            falsePositives = 100*[sum(sortResults(1).invented) ;sum(sortResults(2).invented)]/trueNrSpikes;
            falseNegatives = 100*[sum(~sortResults(1).found) ;sum(~sortResults(2).found)]/trueNrSpikes;
            fprintf('False Positive. Sham: %.0f%%  tACS: %.0f%% \n',falsePositives(1),falsePositives(2))
            fprintf('False Negative. Sham: %.0f%%  tACS: %.0f%% \n',falseNegatives(1),falseNegatives(2))
            fprintf('PLV. Sham: %.2f  tACS: %.2f \n',sqrt(max(0,sortResults(1).ppc)),sqrt(max(0,sortResults(2).ppc)))
            fprintf('PLV2. Sham: %.2f  tACS: %.2f \n',sqrt(max(0,sortResults(1).ppc2)),sqrt(max(0,sortResults(2).ppc2)))

        end
