% grab cell # in its tuning curve across days
% use output from Ele's masterhrz script (modified by Zahra)
% REMEMBER THAT THESE ARE THE ONLY 107 CELLS TRACKED WITH CELLREG ACROSS
% WEEKS AND THEN MAPPED TO DAYS WITH CELLREG

clear all;
fls = dir('Z:\HRZ_master_output8trials\*.mat');
sessions=length(fls);

%%
clear TC_imagesc_arr comp_arr day_arr
cellno=89; %cell to test
for day=1:sessions % grab across days
    load(fullfile(fls(day).folder, fls(day).name))
    fl=fullfile(fls(day).folder, fls(day).name);
    day_cd=fl(end-28:end-25); % hack lol
    comparisons = size(cs_table,1);
    try
    for this_comparison=1:comparisons
        % imagesc plot
        %sorts by the maximum firing rate of the each cell
        [~,max_bin1] = max(cs_table.tuning_curves_for_this_comparison{this_comparison,1},[],2);
        [~,max_bin2] = max(cs_table.tuning_curves_for_this_comparison{this_comparison,2},[],2);
        [~,sorted_idx] = sort(max_bin1);
        
        TC_imagesc = [cs_table.tuning_curves_for_this_comparison{this_comparison,1}(sorted_idx,:) ...
            cs_table.tuning_curves_for_this_comparison{this_comparison,2}(sorted_idx,:)];
        
        % want to collect this across epochs and days...
        TC_imagesc=TC_imagesc(sorted_idx==cellno,:);
    
        figure;
        imagesc(normalize(TC_imagesc,2,'range'))%imagesc(normalize(TC_imagesc,2,'range'))
        colormap(turbo)
        xlabel('cm')
        xticks([0 max([max_bin1; max_bin2])])
        xticklabels([0 180])
        yticks([1])
        yticklabels(cellno) %label with cell id number
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','Times','fontsize',15)
        title(['Day=' day_cd '; TC last 8 trials'],{['epoch ' ...
            num2str(cs_table.comparison_type(this_comparison,1)) ';epoch ' num2str(cs_table.comparison_type(this_comparison,2))]})                      
        if ~sum(TC_imagesc)==0
            TC_imagesc_arr(day,this_comparison,:) = TC_imagesc;
            comp_arr(day,this_comparison,:) = [cs_table.comparison_type(this_comparison,1) cs_table.comparison_type(this_comparison,2)];
            day_arr{day}=day_cd;
        end
    end
    end
end

%%
%test
plt=squeeze(TC_imagesc_arr(:,1,:)); %first comparison
figure;imagesc(normalize(plt,2,'range'));yticks(1:length(day_arr));yticklabels(day_arr)