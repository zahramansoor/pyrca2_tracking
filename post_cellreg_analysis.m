% Zahra
% get cells detected in cellreg and do analysis

clear all
load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230118\Results\cellRegistered_20230118_173248.mat')
% find cells in all sessions
[r,c] = find(cell_registered_struct.cell_to_index_map~=0);
[counts, bins] = hist(r,1:size(r,1));
cindex = bins(counts==3); % finding cells across all 4 sessions yields 0
sessions=4;% specify no of sessions
commoncells=zeros(length(cindex),sessions);
for ci=1:length(cindex)
    commoncells(ci,:)=cell_registered_struct.cell_to_index_map(cindex(ci),:);
end

% load Fall suite2p file
d06=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230118\221206_Fall.mat');
d07=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230118\221207_Fall.mat');
d08=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230118\221208_Fall.mat');
d09=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230118\221209_Fall.mat');
days = [d06,d07,d08,d09];
% calculate dff across all
for i=1:sessions
    dff{i} = redo_dFF(days(i).F, 31.25, 20, days(i).Fneu);
    disp(i)
end

%%% test commoncell#1 - plot across sessions
% plot all cells with legend
close all; 
for cll=1:length(commoncells) % for all cells
    disp(cll)
    figure;
    for i=1:sessions
        if commoncells(cll,i)~=0 % if index is 0 (cell not found in that session)            
            plot(dff{i}(commoncells(cll,i),:));
        else
            yline(2);
        end
        hold on;
    end
    legend('day 1','day 2', 'day 3', 'day 4')
    xlabel('Frames') 
    ylabel('\Delta F / F') 
    savefig(sprintf('C:\\Users\\Han\\Documents\\Zahra\\cellregtest_dff_cell_no_%03d_day_%03d.fig', cll, i))
end
