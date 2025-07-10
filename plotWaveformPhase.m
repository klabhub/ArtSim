function plotWaveformPhase(binCenters,meanWf,sdWf,pv)
% Plot spike waveforms as a function of a phase
arguments
    binCenters double
    meanWf double
    sdWf double
    pv.radius (1,1) double = 0.2
    pv.width(1,1) double = 0.2
    pv.height (1,1) double = 0.2
    pv.alpha (1,1) double = 0.25

    pv.useBins (1,:) double = 1:numel(binCenters)
    pv.flipShading (1,1) logical = false
    pv.colors (1,2) char = 'rg'
    pv.ax = []
    pv.legendBin = []
    pv.legend = {}
    pv.legendLocation (1,1) string = "best"
end


nrPhaseBins = numel(binCenters);
phaseBins = intersect(pv.useBins,1:nrPhaseBins);
t = (1:size(meanWf,1))';

if isempty(pv.ax)
    fig = gcf;
else
    fig = pv.ax.Parent;

end
oldUnits = fig.Units;
fig.Units = 'centimeter';

if ~isempty(pv.ax)

    W =pv.ax.Position(3);
    H = pv.ax.Position(4);
    width  = W.*pv.width;
    height = H.*pv.height;
    radius  = sqrt(sum([W H].^2))*pv.radius;
    center = [pv.ax.Position(1)+0.5*W  pv.ax.Position(2)+0.5*H] -0.5*[width height];

else

    radius = pv.radius;
    width  = pv.width;
    height= pv.height;
    center = [0.5 0.5]-0.5*[width height];
end

pos = [center width height];
h = axes('Position',pos);
for e=phaseBins
    phi = binCenters(e);
    quiver(0,0,cos(phi),sin(phi),'k');
    text(1.2*cos(phi),1.2*sin(phi),num2str(phi*180/pi,3),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',8);
    hold on
end
h.Visible='off';
axs = [];
for e=phaseBins
    phi = binCenters(e);
    pos = [center+radius*[cos(phi) sin(phi)] width height];
    axs = [axs axes('Position',pos)];
    for i=1:2
        if pv.flipShading
            iForErr = setxor(i,[1 2]);
        else
            iForErr = i;
        end
        lower = meanWf(:,e,i)-sdWf(:,e,iForErr);
        upper =  meanWf(:,e,i)+sdWf(:,e,iForErr);
        f= fill([t;t(end:-1:1)], [upper;lower(end:-1:1)],pv.colors(i));
        f.FaceAlpha = pv.alpha;
        f.EdgeAlpha  = 0;
        hold on
        h(i) = plot(t,meanWf(:,e,i),'Color', pv.colors(i),'lineWidth',1);
        axis off
    end
    if ismember(e,pv.legendBin)
        legend(h,pv.legend{:},'Location',pv.legendLocation);
    end
end
linkaxes(axs);
fig.Units = oldUnits;
end

