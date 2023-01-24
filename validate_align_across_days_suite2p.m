% by Zahra
% plot cell mask and activity across days after aligning cells across days in suite2p

% 1- get mean image of each session
% 2- plot cell mask (20/1080 cells) on each of the mean images

d06=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230119\221206_Fall.mat');
d07=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230119\221207_Fall.mat');
d08=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230119\221208_Fall.mat');
d09=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230119\221209_Fall.mat');
d11=load('C:\Users\Han\Documents\MATLAB\CellReg\ZD_Test_20230119\221211_Fall.mat');
days = [d06,d07,d08,d09,d11];

% colormap to iterate thru
ctab = hsv(20);
% randomize cell #
cl=randi(1080,1,20);
figure;
%day 6
subplot(2,3,1) % 2 rows, 3 column, 1 pos
imagesc(d06.ops.meanImg)%imagesc(mat.ops.meanImg);
colormap('gray')
hold on;
for i=1:20
    plot(stat{1,cl(i)}.xpix, stat{1,cl(i)}.ypix, 'Color', ctab(i,:)) % i here is the day
end
axis off
title('day 6')
%day 6
subplot(2,3,2) % 2 rows, 3 column, 1 pos
imagesc(d07.ops.meanImg)%imagesc(mat.ops.meanImg);
colormap('gray')
hold on;
for i=1:20
    plot(stat{1,cl(i)}.xpix, stat{1,cl(i)}.ypix, 'Color', ctab(i,:))
end
axis off
title('day 7')
%day 6
subplot(2,3,3) % 2 rows, 3 column, 1 pos
imagesc(d08.ops.meanImg)%imagesc(mat.ops.meanImg);
colormap('gray')
hold on;
for i=1:20
    plot(stat{1,cl(i)}.xpix, stat{1,cl(i)}.ypix, 'Color', ctab(i,:))
end
axis off
title('day 8')
%day 6
subplot(2,3,4) % 2 rows, 3 column, 1 pos
imagesc(d09.ops.meanImg)%imagesc(mat.ops.meanImg);
colormap('gray')
hold on;
for i=1:20
    plot(stat{1,cl(i)}.xpix, stat{1,cl(i)}.ypix, 'Color', ctab(i,:))
end
axis off
title('day 9')
%day 6
subplot(2,3,5) % 2 rows, 3 column, 1 pos
imagesc(d11.ops.meanImg)%imagesc(mat.ops.meanImg);
colormap('gray')
hold on;
for i=1:20% to account for cell not detected in this day
    plot(stat{1,cl(i)}.xpix, stat{1,cl(i)}.ypix, 'Color', ctab(i,:))
end
axis off
title('day 11')
% savefig(sprintf("C:\\Users\\Han\\Documents\\Zahra\\suite2p_align\\cell_%03d.fig",cl))

savefig(sprintf("C:\\Users\\Han\\Documents\\Zahra\\suite2p_align\\cells.fig"))

%%
% split F.mat across days to fit into cellreg?
bigF = F;
bigFneu = Fneu;
bigspks = spks;
F = bigF(:, 1:45000);
Fneu = bigFneu(:, 1:45000);
spks = bigspks(:, 1:45000);
save('C:\Users\Han\Documents\Zahra\suite2p_align\221206_Fall_mc_across_days.mat', 'F',"Fneu","spks",'ops','iscell','redcell','stat')
F = bigF(:, 45000:90000-1);
Fneu = bigFneu(:, 45000:90000-1);
spks = bigspks(:, 45000:90000-1);
save('C:\Users\Han\Documents\Zahra\suite2p_align\221207_Fall_mc_across_days.mat', 'F',"Fneu","spks",'ops','iscell','redcell','stat')
F = bigF(:, 90000:135000-1);
Fneu = bigFneu(:, 90000:135000-1);
spks = bigspks(:, 90000:135000-1);
save('C:\Users\Han\Documents\Zahra\suite2p_align\221208_Fall_mc_across_days.mat', 'F',"Fneu","spks",'ops','iscell','redcell','stat')
F = bigF(:, 135000:180000-1);
Fneu = bigFneu(:, 135000:180000-1);
spks = bigspks(:, 135000:180000-1);
save('C:\Users\Han\Documents\Zahra\suite2p_align\221209_Fall_mc_across_days.mat', 'F',"Fneu","spks",'ops','iscell','redcell','stat')
F = bigF(:, 180000:225000-1);
Fneu = bigFneu(:, 180000:225000-1);
spks = bigspks(:, 180000:225000-1);
save('C:\Users\Han\Documents\Zahra\suite2p_align\221211_Fall_mc_across_days.mat', 'F',"Fneu","spks",'ops','iscell','redcell','stat')