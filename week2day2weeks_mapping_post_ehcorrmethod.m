% Zahra
% get cells detected in cellreg and do analysis

clear all
load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_week2daymap_week1_20230206\Results\cellRegistered_20230206_180900.mat')
% find cells in all sessions
[r,c] = find(cell_registered_struct.cell_to_index_map~=0);
[counts, bins] = hist(r,1:size(r,1));
sessions=5;% specify no of sessions
cindex = bins(counts==sessions); % finding cells across all 5 sessions
commoncells=zeros(length(cindex),sessions);
for ci=1:length(cindex)
    commoncells(ci,:)=cell_registered_struct.cell_to_index_map(cindex(ci),:);
end
 
% load mats from all days
fls = dir('Z:\week2day_mapping_cellreg\week1\*YC_Fall.mat');%dir('Z:\cellreg1month_Fmats\*YC_Fall.mat');
days = cell(1, length(fls));
for fl=1:length(fls)
    day = fls(fl);
    days{fl} = load(fullfile(day.folder,day.name));
end

% find map from week 1 to those tracked across weeks in previous run
% cells found thru cellreg
weeks=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_across4weeks_20230201\Results\cellRegistered_20230201_154819.mat');
% find cells in all sessions
[r,c] = find(weeks.cell_registered_struct.cell_to_index_map~=0);
[counts, bins] = hist(r,1:size(r,1));
sessions=4;% specify no of sessions
cindex = bins(counts==sessions); % finding cells across all 5 sessions
commoncells_4weeks=zeros(length(cindex),sessions);
for ci=1:length(cindex)
    commoncells_4weeks(ci,:)=weeks.cell_registered_struct.cell_to_index_map(cindex(ci),:);
end

% concat week to week map
% in week map, week 1 is the last element (e.g. 5), whereas in the concat
% week map, it is the first element (bc 'week 1')
week1_across_weeks = commoncells_4weeks(:,1);
% only cells that have both a week to day and week to week mapping
for i=1:5
    week1_mapped2days_across_weeks(:,i) = commoncells(ismember(commoncells(:,5), week1_across_weeks),i);
end
save('Z:\week2day_mapping_cellreg\week1\commoncells_that_are_common_across_weeks.mat','week1_mapped2days_across_weeks')
% now these are cells that are from week 1 that are mapped to the cells
% that are common across weeks
% next step is to do this for each week and save individual mats of this
% map; remember that the last element is always the ROIs of the
% concatenated day

%%
% repeat for week 2
load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_week2daymap_week2_20230206\Results\cellRegistered_20230206_183305.mat')
% find cells in all sessions
[r,c] = find(cell_registered_struct.cell_to_index_map~=0);
[counts, bins] = hist(r,1:size(r,1));
sessions=8;% specify no of sessions
cindex = bins(counts==sessions); % finding cells across all 5 sessions
commoncells=zeros(length(cindex),sessions);
for ci=1:length(cindex)
    commoncells(ci,:)=cell_registered_struct.cell_to_index_map(cindex(ci),:);
end
 
% load mats from all days
fls = dir('Z:\week2day_mapping_cellreg\week2\*YC_Fall.mat');%dir('Z:\cellreg1month_Fmats\*YC_Fall.mat');
days = cell(1, length(fls));
for fl=1:length(fls)
    day = fls(fl);
    days{fl} = load(fullfile(day.folder,day.name));
end

week2_across_weeks = commoncells_4weeks(:,2);
% only cells that have both a week to day and week to week mapping
for i=1:sessions
    week2_mapped2days_across_weeks(:,i) = commoncells(ismember(commoncells(:,sessions), week2_across_weeks),i);
end
save('Z:\week2day_mapping_cellreg\week2\commoncells_that_are_common_across_weeks.mat','week2_mapped2days_across_weeks')

% repeat for week 3
load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_week2daymap_week3_20230203\Results\cellRegistered_20230206_150749.mat')
% find cells in all sessions
[r,c] = find(cell_registered_struct.cell_to_index_map~=0);
[counts, bins] = hist(r,1:size(r,1));
sessions=4;% specify no of sessions
cindex = bins(counts==sessions); % finding cells across all 5 sessions
commoncells=zeros(length(cindex),sessions);
for ci=1:length(cindex)
    commoncells(ci,:)=cell_registered_struct.cell_to_index_map(cindex(ci),:);
end
 
% load mats from all days
fls = dir('Z:\week2day_mapping_cellreg\week3\*YC_Fall.mat');%dir('Z:\cellreg1month_Fmats\*YC_Fall.mat');
days = cell(1, length(fls));
for fl=1:length(fls)
    day = fls(fl);
    days{fl} = load(fullfile(day.folder,day.name));
end

week3_across_weeks = commoncells_4weeks(:,3);
% only cells that have both a week to day and week to week mapping
for i=1:sessions
    week3_mapped2days_across_weeks(:,i) = commoncells(ismember(commoncells(:,sessions), week3_across_weeks),i);
end
save('Z:\week2day_mapping_cellreg\week3\commoncells_that_are_common_across_weeks.mat', ...
    'week3_mapped2days_across_weeks')

% repeat for week 4
load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_week2daymap_week4_20230203\Results\cellRegistered_20230206_143955.mat')
% find cells in all sessions
[r,c] = find(cell_registered_struct.cell_to_index_map~=0);
[counts, bins] = hist(r,1:size(r,1));
sessions=7;% specify no of sessions
cindex = bins(counts==sessions); % finding cells across all 5 sessions
commoncells=zeros(length(cindex),sessions);
for ci=1:length(cindex)
    commoncells(ci,:)=cell_registered_struct.cell_to_index_map(cindex(ci),:);
end
 
% load mats from all days
fls = dir('Z:\week2day_mapping_cellreg\week4\*YC_Fall.mat');%dir('Z:\cellreg1month_Fmats\*YC_Fall.mat');
days = cell(1, length(fls));
for fl=1:length(fls)
    day = fls(fl);
    days{fl} = load(fullfile(day.folder,day.name));
end

week4_across_weeks = commoncells_4weeks(:,4);
% only cells that have both a week to day and week to week mapping
for i=1:sessions
    week4_mapped2days_across_weeks(:,i) = commoncells(ismember(commoncells(:,sessions), week4_across_weeks),i);
end
save('Z:\week2day_mapping_cellreg\week4\commoncells_that_are_common_across_weeks.mat', ...
    'week4_mapped2days_across_weeks')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
week1=load('Z:\week2day_mapping_cellreg\week1\commoncells_that_are_common_across_weeks.mat').week1_mapped2days_across_weeks;
week2=load('Z:\week2day_mapping_cellreg\week2\commoncells_that_are_common_across_weeks.mat').week2_mapped2days_across_weeks;
week3=load('Z:\week2day_mapping_cellreg\week3\commoncells_that_are_common_across_weeks.mat').week3_mapped2days_across_weeks;
week4=load('Z:\week2day_mapping_cellreg\week4\commoncells_that_are_common_across_weeks.mat').week4_mapped2days_across_weeks;
weeks=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_across4weeks_20230201\Results\cellRegistered_20230201_154819.mat');

% find cells in all sessions
[r,c] = find(weeks.cell_registered_struct.cell_to_index_map~=0);
[counts, bins] = hist(r,1:size(r,1));
sessions=4;% specify no of sessions
cindex = bins(counts==sessions); % finding cells across all 5 sessions
commoncells_4weeks=zeros(length(cindex),sessions);
for ci=1:length(cindex)
    commoncells_4weeks(ci,:)=weeks.cell_registered_struct.cell_to_index_map(cindex(ci),:);
end

% for each week, we need to find cells that 1) map to days of that week
% (already in week1,2,... arrays) and 2) map to other days of the week
% need logicals i.e cell 1 in week 1 is in week 2 and week 3 and week 4
% see below...
week1cells_to_map=commoncells_4weeks(:,1); % start with all cells across weeks
sessions_total=20;
cellmap2dayacrossweeks=zeros(length(week1cells_to_map),sessions_total);
for w=1:length(week1cells_to_map)
    week1cell=week1cells_to_map(w);
    daysweek1cell=week1(find(week1(:,end)==week1cell),1:end-1);
    cell_across_weeks=commoncells_4weeks(find(commoncells_4weeks(:,1)==week1cell),:);
    daysweek2cell=week2(find(week2(:,end)==cell_across_weeks(2)),1:end-1); % 1:end-1 to remove week column
    daysweek3cell=week3(find(week3(:,end)==cell_across_weeks(3)),1:end-1);
    daysweek4cell=week4(find(week4(:,end)==cell_across_weeks(4)),1:end-1);
    if ~isempty(daysweek1cell) && ~isempty(daysweek2cell) && ~isempty(daysweek3cell) && ~isempty(daysweek4cell) %make sure cell exists across all days
        cellmap2dayacrossweeks(w,:)=[daysweek1cell,daysweek2cell,daysweek3cell,daysweek4cell];
    end
end

% figures for validation
% align each common cells across all days with an individual mask
% remember this is the cell index, so you have to find the cell in the
% original F mat
cc=cellmap2dayacrossweeks(all(cellmap2dayacrossweeks,2),:);%only gets non zero elements
ctab = hsv(length(cc));

% load mats from all days
fls = dir('Z:\cellreg1month_Fmats\*YC_Fall.mat');%dir('Z:\cellreg1month_Fmats\*YC_Fall.mat');
days = cell(1, length(fls));
for fl=1:length(fls)
    day = fls(fl);
    days{fl} = load(fullfile(day.folder,day.name));
end

for i=1:length(cc)
    %multi plot of cell mask across all 5 days
    figure(i); 
    axes=cell(1,sessions_total);
    for ss=1:sessions_total        
        day=days(ss);day=day{1};
        axes{ss}=subplot(4,5,ss); % 2 rows, 3 column, 1 pos; 20 days
        imagesc(day.ops.meanImg) %meanImg or max_proj
        colormap('gray')
        hold on;
        try
            %plot(day.stat{1,cc(i,ss)}.xpix, day.stat{1,cc(i,ss)}.ypix, 'Color', [ctab(i,:) 0.3]);
        end
        axis off
        title(sprintf('day %i', ss)) %sprintf('day %i', ss)
        %title(axes{ss},sprintf('Cell %0d4', i))
    end
    linkaxes([axes{:}], 'xy')
    %savefig(sprintf("Z:\\suite2pconcat1month_commoncellmasks\\cell_%03d.fig",i+250)) %changed to reflect subset of cells plotted
end

%%
% align all cells across all days in 1 fig
% colormap to iterate thru
ctab = hsv(length(cc));
figure;
axes=zeros(1,sessions_total);
for ss=1:sessions_total
    day=days(ss);day=day{1};
    axes(ss)=subplot(4,5,ss);%(4,5,ss); % 2 rows, 3 column, 1 pos; 20 days
    imagesc(day.ops.meanImg)
    colormap('gray')
    hold on;
    for i=1:length(commoncells)
        try
            plot(day.stat{1,cc(i,ss)}.xpix, day.stat{1,cc(i,ss)}.ypix, 'Color', [ctab(i,:) 0.3]);
        end
    end
    axis off
    title(sprintf('day %i', ss))
end
linkaxes(axes, 'xy')
%savefig(sprintf("Z:\\202300201cells.fig"))

load('Z:\\dff_221206-30.mat')

%%
% plot F (and ideally dff) over ypos

days_to_plot=[1,6,10,14,18]; %plot 5 days at a time
cellno=23;
grayColor = [.7 .7 .7];
fig=figure;
ax1=subplot(5,1,1);
day=days(days_to_plot(1));day=day{1};
plot(day.ybinned, 'Color', grayColor); hold on; 
plot(day.changeRewLoc, 'b')
plot(find(day.licks),day.ybinned(find(day.licks)),'r.')
yyaxis right
try
    plot(dff{days_to_plot(1)}(cc(cellno,days_to_plot(1)),:),'g') % 2 in the first position is cell no
end
title(sprintf('day %i', days_to_plot(1)))

ax2=subplot(5,1,2);
day=days(days_to_plot(2));day=day{1};
plot(day.ybinned, 'Color', grayColor); hold on; 
plot(day.changeRewLoc, 'b')
plot(find(day.licks),day.ybinned(find(day.licks)),'r.')
yyaxis right
try
    plot(dff{days_to_plot(2)}(cc(cellno,days_to_plot(2)),:),'g') % 2 in the second position here is session 2 (day 7)
end
title(sprintf('day %i', days_to_plot(2)))

ax3=subplot(5,1,3);
day=days(days_to_plot(3));day=day{1};
plot(day.ybinned, 'Color', grayColor); hold on; 
plot(day.changeRewLoc, 'b')
plot(find(day.licks),day.ybinned(find(day.licks)),'r.')
yyaxis right
try
    plot(dff{days_to_plot(3)}(cc(cellno,days_to_plot(3)),:),'g') 
end
title(sprintf('day %i', days_to_plot(3)))

ax4=subplot(5,1,4);
day=days(days_to_plot(4));day=day{1};
plot(day.ybinned, 'Color', grayColor); hold on; 
plot(day.changeRewLoc, 'b')
plot(find(day.licks),day.ybinned(find(day.licks)),'r.')
yyaxis right
try
    plot(dff{days_to_plot(4)}(cc(cellno,days_to_plot(4)),:),'g') 
end
title(sprintf('day %i', days_to_plot(4)))

ax5=subplot(5,1,5);
day=days(days_to_plot(5));day=day{1};
plot(day.ybinned, 'Color', grayColor); hold on; 
plot(day.changeRewLoc, 'b')
plot(find(day.licks),day.ybinned(find(day.licks)),'r.')
yyaxis right
try
    plot(dff{days_to_plot(5)}(cc(cellno,days_to_plot(5)),:),'g') 
end
title(sprintf('day %i', days_to_plot(5)))

linkaxes([ax1 ax2 ax3 ax4 ax5],'xy')
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Y position');
xlabel(han,'Frames');
title(han,sprintf('Cell no. %03d', cellno));
savefig(sprintf('Z:\\cellregtest_behavior\\cell_%05d_days%02d_%02d_%02d_%02d_%02d.fig', cellno, days_to_plot))

%%
% align to behavior (rewards and solenoid) for each cell?
% per day, get this data...
range=5;
bin=0.2;
ccbinnedPerireward=cell(1,length(days));
ccrewdFF=cell(1,length(days));
for d=1:length(days)
    day=days(d);day=day{1};
    rewardsonly=day.rewards==1;
    cs=day.rewards==0.5;
    % runs for all cells
    [binnedPerireward,allbins,rewdFF] = perirewardbinnedactivity(dff{d}',rewardsonly,day.timedFF,range,bin); %rewardsonly if mapping to reward
    % now extract ids only of the common cells
    ccbinnedPerireward{d}=binnedPerireward(cc(:,d),:);
    ccrewdFF{d}=rewdFF(:,cc(:,d),:);
end
%%
% plot
%cellno=34; % cell to align
optodays=[5,6,7,9,10,11,13,14,16,17,18];
for cellno=1:length(cc) %or align all cells hehe
    figure;
    for d=1:length(days)
        pltrew=ccbinnedPerireward{d};
        if ~any(optodays(:)==d)            
            plot(pltrew(cellno,:)', 'Color', 'black')            
        else
            plot(pltrew(cellno,:)', 'Color', 'red')    
        end
        hold on;        
    end
    % plot reward location as line
    xticks([1:5:50, 50])
    x1=xline(median([1:5:50, 50]),'-.b','Reward'); %{'Conditioned', 'stimulus'}
    xticklabels([allbins(1:5:end) range]);
    xlabel('seconds')
    ylabel('dF/F')
    title(sprintf('Cell no. %04d', cellno))
end