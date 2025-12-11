function [f,ax] = fig(varargin)
% Generate a figure matching publisher requirements.

defaultColorOrder = [   0.6350    0.0780    0.1840
    0.4660    0.6740    0.1880
    0.9290    0.6940    0.1250
    0    0.4470    0.7410
    0.3010    0.7450    0.9330
    0.4940    0.1840    0.5560
    0.8500    0.3250    0.0980];

p=inputParser;
p.addParameter('paperCols',1,@(x) (ismember(x,[1 1.5 2])));
p.addParameter('byColumn',false,@islogical);
p.addParameter('nrCols',1);
p.addParameter('nrRows',1);
p.addParameter('height',NaN);
p.addParameter('name','');
p.addParameter('colorOrder',defaultColorOrder);
p.parse(varargin{:});


%% Set Defaults
set(groot,'defaultAxesFontSize',8)
set(groot,'defaultAxesFontName','Arial');
set(groot,'defaultAxesLineWidth',0.5);


%% Find or create the figure
f = findall(0,'type','figure','Name',p.Results.name);
if ~isempty(f)
    close(f)
end

f =figure;



switch (p.Results.paperCols)
    case 1
        width = 6.68;
    case 1.5
        width = 9.94;
    case 2
        width = 13.2;
end
if isnan(p.Results.height)
    height = width;
else
    height = p.Results.height;
end
clf(f,"reset");
f.Name = p.Results.name;
f.Units= 'centimeters';
f.Position = [f.Position(1:2) width height];
f.PaperPosition= [f.PaperPosition(1:2) width height];

figure(f);
%% Create the axes
nrSubs = p.Results.nrRows*p.Results.nrCols;
layout = tiledlayout(p.Results.nrRows,p.Results.nrCols); %#ok<NASGU>
for i =1:nrSubs
    ax(i)=nexttile;
    cla(ax(i));
    set(ax(i),'ColorOrder',p.Results.colorOrder);
    hold on;
end
ax = reshape(ax,[p.Results.nrRows p.Results.nrCols]);
if p.Results.byColumn
    ax= ax';
end

for i=1:nrSubs
    if nrSubs>1
        thisPos = ax(i).Position;
        thisPos(1) = max(0,thisPos(1) - 0.35*thisPos(3));
        thisPos(2) = thisPos(2) + 0.95*thisPos(4);
        thisPos(3) = 0.1;
        thisPos(4) = 0.1;
        annotation('textbox',thisPos,'String',[char('a'+i-1) ')'],'Linestyle','none','FontName',get(groot,'defaultAxesFontName'))
    end
end