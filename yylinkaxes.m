function yylinkaxes(axs,dimension,which)
% Linkaxes function that can handle axes that use two y-axes (from yyaxis)
%
arguments
    axs  (1,:)             % Handles to axes that should be linked
    dimension (1,1)  string = "xyz"   % Dimension to be linked 
    which (1,1) string {mustBeMember(which,["both","left","right"])} ="both"  % Yaxis to be linked ("left" "right" or "both")
end

yaxes= {axs.YAxis};
isyy = cellfun(@(x) numel(x)>1,yaxes);

if ~any(isyy)
    linkaxes(axs,dimension);
else
    if which =="both"
        for a=axs
            yyaxis(a,"left")
        end
        linkaxes(axs,dimension);
        for a=axs
            yyaxis(a,"right")
        end
    else
        for a=axs
            yyaxis(a,which)
        end
        linkaxes(axs,dimension);
    end
end