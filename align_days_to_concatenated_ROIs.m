% ***need to fix all_F output, if any 0s, set entire row to NaN


% For concatenating ROIs across days
% 210702 Notes on data pipeline
% run "runVideosTiff_EH_new_sbx_setCrop_planes" to crop FOVs (important for limiting size) and split planes 
% I crop aggressively for each plane. if we lose some cells at edges, that's ok
% calls "loadVideoTiffNoSplit_EH2_new_sbx_setCrop_planes"
% then this file, which also allows you to check alignments

% Aligns cell ROIs across days based on aligning single suite2p day to all days concatenated suite2p data
% Uses individual day ROI fits which should have better S/N than across day ROIs
% 1st, make tifs for each day and plane, then into folders "d1" to "d(x)"
% Within each day folder split planes to folders "d1\p1" to "d1\p4"
% This is important because you don't want to run suite2p on multiple planes at once, just 1!
% Run suite2p on each day and plane individually. 
% not sure we need to save bin or run two step registration for these. however, not tested
% The original folder and subfolders, "d1\p1", etc for day 1, plane1
% have the single day/single plane analyses.
% I do minimal curation on these, just a quick check over.
% Do this for all days and planes. Can do batch mode.
% 
% Make a concatenated suite2p analysis. only one plane at a time. 
% for ex, make "d1-d11_p1" folder.
% Then run suite2p by adding all the individual day, p1 folders to the path
% set the save path to the d1-d11_p1 folder.
% This registers all the tifs from the individual folders together
% This is used as a template to align individual days. the actual dF is not used for data.
% Spend time to check that your parameters are optimal for motion corr and aligning movies.
% two step registration helps and save bin etc to look at motion corr
% In particular "maxregshiftNR" may need to be changed. 
% on CA1 movies across 11 days, i found 10 pix to be optimal. assume this will vary
% as always it's good to run on a shorter test data set to optimize parameters
% other default settings, threshold 0.9, max_iterations 30. gets lots of rois, not all good
% after running, move "suite2p" folder of concatenated movie (will be in "d1\p1" to e.g. "d1-11_p1"
% Check quality and curate these ROIs carefully
% should have individual days and curated concatenated suite2p for all planes
% 
% this script asks you number of days to align and the locations of concatenated and individual day files
% it will save this data in a file list in case you want to re-run the analysis

%then goes through and selects best fits to template based on F, plots matrix and r vs. cent_dist.
% can look at template vs single day ROIs and traces, or skip
% if centroid dist is too far, writes as NaN
%saves file list and file
%moved code to assess ROI fits to separate script that can be run later
%"template_to_single_day_ROI_comparison", however not working


% Saves data in "Fall_xd_align.mat" in suite2p folder of concatenated movie
% file list also save here

close all
clear all

%variables for 
max_cent_dist=20;%maximum centroid dist in pixels, should convert to um based on mag
% in test experiment, was 2X maxshift for non rigid registration (10pix)
min_r=0.2; %min r value of matching roi. probably doesn't matter much
% in testing, centroid distance was the most critical variable for finding borderline matches
% at low r values and distances, may or may not be same cell. 
% if it's the same cell, very little signal anyway

answer = questdlg('Do you have a saved file list?',	'Loading F files and paths','Yes','No','No');
switch answer
    case 'Yes'
        [pre_list,pre_path]=uigetfile('*_file_list.mat','Select file list');
        cd (pre_path); %set path
        load(pre_list);%in concat day folder
    case'No'
        prompt = {'Number of Files:'};
        dlgtitle = 'Days to align';
        dims = [1 35];
        definput = {'5'}; %default files
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        numDays=str2num(answer{1});
        single_files=cell(numDays,1);
        single_paths=cell(numDays,1);

        [all_file,all_path]=uigetfile('*.mat','All days concatenated file ');
        cd (all_path); %set path

        for j=1:numDays
            [single_files{j,1},single_paths{j,1}]=uigetfile('*.mat',['Single day analysis, Day ',num2str(j)]);
            cd (single_paths{j,1}); %set path
        end
end
cd (all_path); %set path
templ_all=load('Fall.mat');%in concat day folder
F_all=templ_all.F';
rois = find(templ_all.iscell(:,1));%Nx1 vector with index of iscell=1, template cells
%preallocate, variables below filled at end of loop
frames=zeros(numDays,1);%keep track of #frames per day
dayStart=zeros(numDays,1);%keep track 1st frame of days
dayEnd=zeros(numDays,1);%keep track of last frame of days
%below filled each loop, prior to filtering for bad fits
all_M=zeros(size(rois,1),numDays);%max r per roi across days
all_I=zeros(size(rois,1),numDays);%index of max r cell for rois across days
all_cent_dist=zeros(size(rois,1),numDays);
%after filtering for bad ROI fits, look up table below
LUT=zeros(size(rois,1),numDays);%look up table of all roi to single across days. not found or poor fit = NaN

%all has stitched single F from cross correlated ROIs to template ROIs
%could do Fc2/3 on individual days or across all. need to think about
%for now, is commented out
all_F= NaN(size(templ_all.F));
all_Fneu= NaN(size(templ_all.F));
all_spks= NaN(size(templ_all.F));
% all_Fc2= NaN(size(templ_all.F));
% all_Fc3= NaN(size(templ_all.F));
%structure with corrcoef results for each all day ROI
% day.corrs(day,1)
%structure with single.stat for each single day, ROIs info
% day.stat(day,1);
endIdx=0;
                                                
% all_F= all_F';%transpose back to original
% all_Fc2= all_Fc2';
% all_Fc3= all_Fc3';

stripped_name=regexprep(all_file,'.mat','');
currfile=strcat(stripped_name,'_%dd_align.mat');
currfile=sprintf(currfile,numDays);
currfilename=[all_path currfile];
save(currfilename,'-v7.3');    %need -v7.3 MAT file or variable is too big to save

file_list=strcat(stripped_name,'_%dd_file_list.mat');
file_list=sprintf(file_list,numDays);
file_list_full=[all_path file_list];
save(file_list_full,'all_file','all_path','single_files','single_paths','numDays', '-v7.3');    %need -v7.3 MAT file or variable is too big to save


