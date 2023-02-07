

for i =1:numel(LUT(:,1))%loop through all rois in column 1
    dupes(i,1)=numel(find(LUT(:,1)==LUT(i,1)));%find total instances in column with that roi
 end

dupIdx=find(dupes(:,1)>1);%idx of duplicates. nan=0, uniques=1, duplicates>1
dupROIs=LUT(dupIdx,1);