function recreate_iscell_and_make_all_struct(Settings)

for this_day = 1:size(Settings.paths,1)

    clearvars -except this_day Settings
%     should_be_analyzed = 1;

    file =fullfile(Settings.paths(this_day).folder,Settings.paths(this_day).name);
%     directory = file;
%     info = split(directory,'/');
%     mouse_cd = string(info{Settings.level_mouse_name});
%     day_cd = string(info{Settings.level_day});

   l = load(file);
   %skewdcells = find(skewness(l.F,1,2)<2); %looks at skewness of cells, <2 --> interneurons; if omitted, gets all cells
   %l.iscell(skewdcells,1) = 0;
   disp ([file ' ... loaded'])
    %don't remove iscells==0, unnecessary for this analysis
    try
        zdtest=l.all; % ZD added
    end
    if exist('zdtest', 'var')
        all=l.all; % ZD added
    end
    if ~exist('all','var')
        disp('problem with all variable ')
        disp('..recreating all structure.. ')
        [x,y]=size(l.F);
        cells = ones(1,x) ; %take all cells

        cd(Settings.paths(this_day).folder)
        all = create_all_structure(l.F,l.Fneu,l.spks,Settings.Fs,cells);
        save(file , 'all','-append') 
        disp('done recreating all !')
    end
    disp ([file ' ... done!'])
end
end