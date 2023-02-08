
% plot from VR file
ypos = VR.ypos;
ypos(ypos<2) = nan;
time_min = VR.time/60;
time_min(ypos<2) = nan;
licks = VR.lick;
rewards = VR.reward;
changeRewLoc = find(VR.changeRewLoc);
RewLoc = VR.changeRewLoc(changeRewLoc>0);
Gain = VR.scalingFACTOR;


figure;
plot(time_min,ypos,'k')
hold on
scatter(time_min(licks==1),ypos(licks==1),15,"filled",'r')
scatter(time_min(rewards==1),ypos(rewards==1),15,"filled",'b')
for mm = 1:length(changeRewLoc)-1 %the rectangle indicating the reward location, overlaps the probe trials referring to the previous reward location
    rectangle('position',[time_min(changeRewLoc(mm)) RewLoc(mm)-7.5*Gain time_min(changeRewLoc(mm+1))-time_min(changeRewLoc(mm)) 15*Gain],'EdgeColor',[0 0 0 0],'FaceColor',[0 .5 .5 0.3])
end
ylabel('track position - cm')
xlabel('time_min - minutes')
title ('lick (red) along track - and reward (blue)')