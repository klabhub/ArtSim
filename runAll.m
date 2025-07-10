function runAll(trgFoldr, exportFormats,mlxFiles)
%% Run all mlx files
% This is intended as a test that all code runs. This will not update the
% MLX files but does generate regular figures, exports them as pnd and pdf
% and the command line output from the mlx files is stored ina a txt file.  
% This output is saved to the trgFolder.
arguments
    trgFoldr    (1,1) string    = "c:/temp/artsim/";
    exportFormats  (1,:) string  = ["-dpng" "-dpdf"];
    mlxFiles      struct = dir('*.mlx');
end

close all;
if ~exist(trgFoldr,"dir")
    mkdir(trgFoldr);
end
for i=1:numel(mlxFiles)
    try
<<<<<<< HEAD
        clearvars -except mlxFiles i trgFoldr exportFormats % Avoid out of memory / scope issues
        name = 
        diary(fullfile(trgFoldr,[extractBefore(mlxFiles(i).name,'.mlx') '.txt']))
=======
        clearvars -except mlxFiles i trgFoldr exportFormats % Avoid out of memory / scope issues        
        name = extractBefore(mlxFiles(i).name,'.mlx');
        fprintf(2,'**************** %s ****************\n',name)
        diary(fullfile(trgFoldr,[name '.txt']))
>>>>>>> 8d8e22a9767cf568b8a6dbc14e7590bd817c9d36
        tic
        run(mlxFiles(i).name)
        figures = findobj(0,'type','figure');
        savefig(fig, fullfile(trgFoldr,extractBefore(mlxFiles(i).name,'.mlx'))); % All figures in one .fig file
        for exportFormat = exportFormats 
            for fig =figures'
                filename = fullfile(trgFoldr,[name get(fig,'Name')]);
                print(fig,exportFormat,filename);
            end
        end
        close(figures)
        clear figures
<<<<<<< HEAD
        
=======
        save(fullfile(trgFoldr,name)) % Save state to .mat file nanmed after mlx
>>>>>>> 8d8e22a9767cf568b8a6dbc14e7590bd817c9d36
        toc
        diary off
    catch me
        fprintf(2,'\n  !!!! %s  !!!! \n',me.message)
        close all
        diary off
    end
end

