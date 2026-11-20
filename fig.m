function [f,ax] = fig(pv)
% Generate a figure matching publisher requirements.

arguments
    pv.paperCols (1,1) double {mustBeMember(pv.paperCols,[1 1.5 2])} = 1
    pv.byColumn (1,1) logical = false
    pv.nrCols (1,1) double = 1
    pv.nrRows (1,1) double = 1
    pv.height (1,1) double = NaN
    pv.name = ''
    pv.colorOrder (:,3) double = [0.6350 0.0780 0.1840
                                   0.4660 0.6740 0.1880
                                   0.9290 0.6940 0.1250
                                   0.0000 0.4470 0.7410
                                   0.3010 0.7450 0.9330
                                   0.4940 0.1840 0.5560
                                   0.8500 0.3250 0.0980]
    pv.span = {}   
    pv.fontSize = 6;
    pv.fontName = 'Arial';
    pv.lineWidth = 1;
    pv.dx = 1; % Number of expected chars on the yticklabels (used to left shift the subplot labels)
end


%% Find or create the figure
f = findall(0,'type','figure','Name',pv.name);
if ~isempty(f)
    close(f)
end

f =figure;

switch (pv.paperCols)
    case 1
        width = 8.8;
    case 1.5
        width = 13.2;
    case 2
        width = 18;
end
if isnan(pv.height)
    height = width;
else
    height = pv.height;
end
clf(f,"reset");
f.Name = pv.name;
f.Units= 'centimeters';
f.Position = [f.Position(1:2) width height];
f.PaperPosition= [f.PaperPosition(1:2) width height];

figure(f);
%% Create the axes
if isempty(pv.span)
    nrSubs = pv.nrRows*pv.nrCols;
else
    nrSubs = numel(pv.span);
end
layout = tiledlayout(pv.nrRows,pv.nrCols); %#ok<NASGU>
ax = gobjects(1,nrSubs);
for i =1:nrSubs
    if isempty(pv.span)
        ax(i)=nexttile;
    else
        ax(i)=nexttile(pv.span{i}(1),pv.span{i}(2:3));
    end
    cla(ax(i));
    set(ax(i), ...
        'ColorOrder',pv.colorOrder, ...
        'FontSize',pv.fontSize, ...
        'LabelFontSizeMultiplier',1, ...
        'TitleFontSizeMultiplier',1);
    hold on;
end

if isempty(pv.span)
ax = reshape(ax,[pv.nrCols pv.nrRows])';
end
% Create labels to the left of each axes, centered at ylabel height.
if nrSubs > 1
    if pv.byColumn
        labelAxes = ax';
    else
        labelAxes = ax;
    end

    for i = 1:nrSubs
        yLabel = labelAxes(i).YLabel;
        yLabel.Units = "normalized";
        hTmp1 = plot(labelAxes(i),0:0.1:10^(pv.dx-1));
        hTmp2 = ylabel(labelAxes(i),"BlaBlaBla");       
        
          % Align the subplot label with the estimated left edge of the ylabel.
        labelX = hTmp2.Position(1) - hTmp2.Extent(3);
        delete([hTmp1 hTmp2]);

      
        text(labelAxes(i), labelX, 1, ...
            string(char('a' + i - 1)) + ")", ...
            Units="normalized", ...
            HorizontalAlignment="left", ...
            VerticalAlignment="bottom", ...
            Clipping="off", ...
            Interpreter="none", ...
            FontName=pv.fontName, ...
            FontSize=pv.fontSize);
    end
end
