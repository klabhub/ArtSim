function runAll(pv)
%% Run all mlx files
% This is intended as a test that all code runs. This will not update the
% MLX files but does generate regular figures, exports them as pnd and pdf
% and the command line output from the mlx files is stored ina a txt file.
% This output is saved to the trgFolder.
arguments
    pv.trgFoldr    (1,1) string    = "c:/temp/artsim/";
    pv.exportFormats  (1,:) string  = ["eps" "pdf"];
    pv.mlxFiles      struct = dir('*.mlx');
    pv.overwrite (1,1) logical = false
end

close all;
if ~exist(pv.trgFoldr,"dir")
    mkdir(pv.trgFoldr);
end
for i=1:numel(pv.mlxFiles)
    try
        clearvars -except pv i
        name = extractBefore(pv.mlxFiles(i).name,'.mlx');
        trgFile = fullfile(pv.trgFoldr,name);
        if ~exist(trgFile + ".mat" ,'file') || pv.overwrite
            fprintf(2,'**************** %s ****************\n',name)
            diary(trgFile + ".txt");
            tic
            run(pv.mlxFiles(i).name)
            figures = findobj(0,'type','figure');
            savefig(figures, trgFile); % All figures in one .fig file
            for exportFormat = pv.exportFormats
                for fg =figures'
                    filename = fullfile(pv.trgFoldr,[name get(fg,'Name')]);
                    exportgraphics(fg,filename + "." + exportFormat,contentType='vector',Colorspace="rgb",PreserveAspectRatio="on",Resolution=600)                    
                end
            end
            close all
            clear figures
            save(trgFile) % Save state to .mat file nanmed after mlx
            toc
            diary off
        else
            fprintf(2,'**************** %s (Skipped) *******\n',name)
        end
    catch me
        fprintf(2,'\n %s: !!!! %s  !!!! \n',name, me.message)
        me.stack
        close all
        diary off
    end
end

