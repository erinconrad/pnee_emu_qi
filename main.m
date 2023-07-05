clear
clc

%% To dos
%{
- Decide what to do with nans
    - e.g., in known_to_have_pnee_prior
- Check that my calculation of the U statistic is correct
%}

%% Dependencies
% This requires adding myFisher23 to your path 
% https://www.mathworks.com/matlabcentral/fileexchange/15399-myfisher23

%% Get file locations and load redcap file
locations = pnee_locations;
data_folder = locations.data;
results_folder = locations.results;
file_name = 'PNEEEMUPathway-Erinanalysisfull_DATA_2023-07-02_0630.csv';
label_file_name = 'PNEEEMUPathway-Erinanalysisfull_DATA_LABELS_2023-07-02_0813.csv';

% load redcap output file into table
T = readtable([data_folder,file_name]);
lT = readtable([data_folder,label_file_name]);

% Get table variable names
var_names_T = T.Properties.VariableNames;
var_names_lT = lT.Properties.VariableNames;

%% Identify pre-intervention and post-intervention patients
all_pts = ones(length(T.phase),1); all_pts = logical(all_pts);
pre = T.phase == 1;
post = T.phase == 2;

%% Group variables according to type
% This will let me know what analysis to run
continuous_normalish = {'age','since_discharge_event_freq',...
    'since_discharge_event_seve','work','home','social','private',...
    'relationships','qol','work1','home1','social1','private1',...
    'relationships1','qol1','since_discharge_event_freq1','since_discharge_event_seve1'};
continuous_skewed = {'event_frequency_pre','follow_up_interval',...
    'event_frequency_post','event_frequency_post1'};
binary = {'known_to_have_pnee_prior','dual_diagnosis',...
    'follows_with_psych_pre','asm_pre','asm_dc','were_asm_for_an_indication',...
    'talked_to_patient','patient_arranged_psych','study_team_member_arranged',...
    'did_patient_follow_up_with','was_follow_up_scheduled',...
    'patient_understands_diagno','did_your_understanding_of',...
    'patient_agrees_with_diagno','ed_visits_or_hospitalizati','patient_scheduled_appointm',...
    'anti_seizure_medication_wa','psychiatry_consult_obtaine','follow_up_with_psychiatris',...
    'psych_med_was_changes_star','self_help_or_apps_used_for','patient_understands_diagno1',...
    'did_your_understanding_of1','patient_agrees_with_diagno1',...
    'ed_visits_or_hospitalizati1',...
    'patient_scheduled_appointm1','anti_seizure_medication_wa1','follow_up_with_psychiatris1',...
    'psych_med_was_changes_star1','self_help_or_apps_used_for1'};
non_binary = {'gender'};

%% Prepare structure that will contain stats for every variable

% Start a running count of all variables
var_count = 0;

%% Fill up continuous normal comparisons in structure
% For now I'll do non parametric stats for consistency, but we can change
% this
ncontinuous_normal = length(continuous_normalish);

% Loop over the continuous normal variables
for i = 1:ncontinuous_normal
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = continuous_normalish{i};
    out.var(var_count).type = 'continuous normalish';

    % Get pre and post for this variable
    var = T.(continuous_normalish{i});
    var_pre = var(pre);
    var_post = var(post);

    % Get median and IQR for pre and post
    pre_median = nanmedian(var_pre);
    pre_iqr = prctile(var_pre,[25 75]);
    post_median = nanmedian(var_post);
    post_iqr = prctile(var_post,[25 75]);

    % Do Wilcoxon rank sum (same as Mann Whitney U) to do non-parametric
    % comparison of two independent groups
    try
        [p,~,stats] = ranksum(var_pre,var_post);
        W = stats.ranksum;
    catch ME
        if contains(ME.message,'No data remaining after removal of NaNs')
            p = nan;
            stats = [];
            W = nan;

        end
    end

    % Get U statistic
    U = W-sum(~isnan(var_pre))*(sum(~isnan(var_pre))+1)/2;

    % Get the column number of this table
    col_num = strcmp(var_names_T,continuous_normalish{i});

    % Get the corresponding variable name from the label table
    lt_variable = var_names_lT{col_num};

    % Prepare text
    stats_text = {sprintf('%s: median (IQR)',lt_variable),...
        sprintf('%1.1f (%1.1f-%1.1f)',pre_median,pre_iqr(1),pre_iqr(2)),...
        sprintf('%1.1f (%1.1f-%1.1f)',post_median,post_iqr(1),post_iqr(2)),...
        sprintf('U: %1.1f',U),...
        sprintf('%s',formatted_p_values(p))};

    % Add this info to the struct
    out.var(var_count).pre_median = pre_median;
    out.var(var_count).pre_iqr = pre_iqr;
    out.var(var_count).post_median = post_median;
    out.var(var_count).post_iqr = post_iqr;
    out.var(var_count).test = 'ranksum';
    out.var(var_count).p = p;
    out.var(var_count).stats = stats;
    out.var(var_count).U = U;
    out.var(var_count).text = stats_text;

end

%% Fill up continuous skewed comparisons in structure
% This should be non parametric comparisons
ncontinuous_skewed = length(continuous_skewed);

% Loop over the continuous normal variables
for i = 1:ncontinuous_skewed
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = continuous_skewed{i};
    out.var(var_count).type = 'continuous skewed';

    % Get pre and post for this variable
    var = T.(continuous_skewed{i});
    var_pre = var(pre);
    var_post = var(post);

    % Get median and IQR for pre and post
    pre_median = nanmedian(var_pre);
    pre_iqr = prctile(var_pre,[25 75]);
    post_median = nanmedian(var_post);
    post_iqr = prctile(var_post,[25 75]);

    % Do Wilcoxon rank sum (same as Mann Whitney U) to do non-parametric
    % comparison of two independent groups
    try
        [p,~,stats] = ranksum(var_pre,var_post);
        W = stats.ranksum;
    catch ME
        if contains(ME.message,'No data remaining after removal of NaNs')
            p = nan;
            stats = [];
            W = nan;
        end
    end

    % Get U statistic. DOUBLE CHECK THAT THIS IS RIGHT
    U = W-sum(~isnan(var_pre))*(sum(~isnan(var_pre))+1)/2;

    % Get the column number of this table
    col_num = strcmp(var_names_T,continuous_skewed{i});

    % Get the corresponding variable name from the label table
    lt_variable = var_names_lT{col_num};

    % Prepare text
    stats_text = {sprintf('%s: median (IQR)',lt_variable),...
        sprintf('%1.1f (%1.1f-%1.1f)',pre_median,pre_iqr(1),pre_iqr(2)),...
        sprintf('%1.1f (%1.1f-%1.1f)',post_median,post_iqr(1),post_iqr(2)),...
        sprintf('U: %1.1f',U),...
        sprintf('%s',formatted_p_values(p))};


    % Add this info to the struct
    out.var(var_count).pre_median = pre_median;
    out.var(var_count).pre_iqr = pre_iqr;
    out.var(var_count).post_median = post_median;
    out.var(var_count).post_iqr = post_iqr;
    out.var(var_count).test = 'ranksum';
    out.var(var_count).p = p;
    out.var(var_count).stats = stats;
    out.var(var_count).text = stats_text;



end

%% Fill up binary comparisons in structure
% This should be Fisher exact tests
nbinary = length(binary);

% Loop over the binary variables
for i = 1:nbinary
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = binary{i};
    out.var(var_count).type = 'binary';

    % Get pre and post for this variable
    var = T.(binary{i});
    var_pre = var(pre);
    var_post = var(post);

    % Get the column number of this table
    col_num = strcmp(var_names_T,binary{i});

    % Get the corresponding column from the label table
    lt_variable = var_names_lT{col_num};
    lT_col = lT.(lt_variable);

    % double check that it really is binary
    assert(sum(var == 1 | var == 0 | isnan(var)) == length(var))

    % get the labels corresponding to 1 and 0
    label1 = unique(lT_col(var==1));
    label0 = unique(lT_col(var==0));
    assert(length(label1)<=1); assert(length(label0)<=1);

    try
        out.var(var_count).labels = {label0{1} label1{1}};
    catch
        if isempty(label0) && strcmp(label1{1},'Yes')
            label0{1} = 'No';
        elseif isempty(label1) && strcmp(label0{1},'No')
            label1{1} = 'Yes';
        end
        out.var(var_count).labels = {label0{1} label1{1}};
    end


    % Get numbers of 1s and 0s for pre and post
    pre1 = sum(var_pre==1);
    pre0 = sum(var_pre==0);
    post1 = sum(var_post == 1);
    post0 = sum(var_post==0);

    % Prep 2x2 table
    tbl_2x2 = table([pre0;pre1],[post0;post1],...
        'VariableNames',{'Pre','Post'},'RowNames',{'0','1'});
    [~,p,stats] = fishertest(tbl_2x2);

    % Prep text
    stats_text = {sprintf('%s',lt_variable),'','',...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label0{1}),...
        sprintf('%d (%1.1f%%)',pre0,pre0/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post0,post0/length(var_post)*100),'','';...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',pre1,pre1/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post1,post1/length(var_post)*100),'',''};

    % Add this to structure
    out.var(var_count).tbl = tbl_2x2;
    out.var(var_count).p = p;
    out.var(var_count).stats = stats;
    out.var(var_count).text = stats_text;

end

%% Fill up non-binary comparisons
nnon_binary = length(non_binary);
% Loop over the non binary variables
for i = 1:nnon_binary
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = non_binary{i};
    out.var(var_count).type = 'non-binary';

    % Get pre and post for this variable
    var = T.(non_binary{i});
    var_pre = var(pre);
    var_post = var(post);

    % double check that options are just 1, 2, and 3
    assert(sum(var <= 1 | var <= 2 | var <= 3) == length(var))

    % Get the column number of this table
    col_num = strcmp(var_names_T,non_binary{i});

    % Get the corresponding column from the label table
    lt_variable = var_names_lT{col_num};
    lT_col = lT.(lt_variable);

    % get the labels corresponding to 1 and 0
    label1 = unique(lT_col(var==1));
    label2 = unique(lT_col(var==2));
    label3 = unique(lT_col(var==3));
    assert(length(label1)==1); assert(length(label2)==1); assert(length(label3)==1);
    out.var(var_count).labels = {label1{1} label2{1} label3{1}};

    % Get numbers of 1s and 0s for pre and post
    pre1 = sum(var_pre==1);
    pre2 = sum(var_pre==2);
    pre3 = sum(var_pre==3);
    post1 = sum(var_post == 1);
    post2 = sum(var_post==2);
    post3 = sum(var_post==3);

    % Prep 2x3 table
    tbl_2x3 = table([pre1;pre2;pre3],[post1;post2;post3],...
        'VariableNames',{'Pre','Post'},'RowNames',{'1','2','3'});
    
    % Convert this into the format needed for MyFisher23
    x = table2array(tbl_2x3);
    x = x';

    % do myfisher23
    p = myfisher23(x);

    % Prep text
    stats_text = {sprintf('%s',lt_variable),'','',...
        '',...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',pre1,pre1/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post1,post1/length(var_post)*100),'','';...
        sprintf('%s',label2{1}),...
        sprintf('%d (%1.1f%%)',pre2,pre2/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post2,post2/length(var_post)*100),'','';...
        sprintf('%s',label3{1}),...
        sprintf('%d (%1.1f%%)',pre3,pre3/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post3,post3/length(var_post)*100),'',''};

    % Add this to structure
    out.var(var_count).tbl = tbl_2x3;
    out.var(var_count).p = p;
    out.var(var_count).text = stats_text;

end

%% Dump all the pre vs post comparisons into one huge table
all_pre_vs_post_comp = [];
n_pre_vs_post_comp = length(out.var);
for i = 1:n_pre_vs_post_comp
    all_pre_vs_post_comp = [all_pre_vs_post_comp;...
        out.var(i).text];
end

%% Turn this into a table and save it as a csv file
outT = cell2table(all_pre_vs_post_comp,'VariableNames',...
    {'Variable','Pre','Post','Statistic','p-value'});
writetable(outT,[results_folder,'pre_vs_post_table.csv']);

%% Prepare paired analyses
pre_freq_var = 'event_frequency_pre';
post_freq_vars = {'event_frequency_post','event_frequency_post1'};
npost = length(post_freq_vars);

% Loop over two ways of defining post-EMU event frequency (visit vs survey)
for i = 1:npost
    
    % Get pre- and post-EMU 
    pre_emu = T.(pre_freq_var);
    post_emu = T.(post_freq_vars{i});
    paired.post_emu_def(i).name = post_freq_vars{i};

    paired.post_emu_def(i).text = {sprintf('Pre-vs-post EMU comparisons using %s to define post-EMU event frequency',...
        post_freq_vars{i}),'','','',''};
    
    % Do for three groups: all patients, just pre-intervention arm, just
    % post-intervention arm
    for isub = 1:3
        if isub == 1
            group = all_pts;
            which_group = 'all patients';
        elseif isub == 2
            group = pre;
            which_group = 'pre-intervention';
        elseif isub == 3
            group = post;
            which_group = 'post-intervention';
        end

        % Get median and IQR for pre and post
        pre_median = nanmedian(pre_emu(group));
        pre_iqr = prctile(pre_emu(group),[25 75]);
        post_median = nanmedian(post_emu(group));
        post_iqr = prctile(post_emu(group),[25 75]);

         % Do sign rank test (paired non-parametric comparison) 
         try
            [p,~,stats] = signrank(pre_emu(group),post_emu(group));
            signedrank = stats.signedrank;
         catch
             p = nan;
             stats = [];
             signedrank = nan;
         end


        % Prepare text
        stats_text = {sprintf('%s',which_group),'','','','';...
            sprintf('Event frequency: median (IQR)'),...
        sprintf('%1.1f (%1.1f-%1.1f)',pre_median,pre_iqr(1),pre_iqr(2)),...
        sprintf('%1.1f (%1.1f-%1.1f)',post_median,post_iqr(1),post_iqr(2)),...
        sprintf('T+: %1.1f',signedrank),...
        sprintf('%s',formatted_p_values(p))};

        % Save stuff
        paired.post_emu_def(i).group(isub).name = which_group;
        paired.post_emu_def(i).group(isub).pre_median = pre_median;
        paired.post_emu_def(i).group(isub).post_median = post_median;
        paired.post_emu_def(i).group(isub).pre_iqr = pre_iqr;
        paired.post_emu_def(i).group(isub).post_iqr = post_iqr;
        paired.post_emu_def(i).group(isub).stats = stats;
        paired.post_emu_def(i).group(isub).p = p;
        paired.post_emu_def(i).group(isub).numbers = [pre_emu(group) post_emu(group)];
        paired.post_emu_def(i).group(isub).text = stats_text;
        if isub == 1
            paired.post_emu_def(i).group(isub).pre_intervention = pre;
            paired.post_emu_def(i).group(isub).post_intervention = post;
        end


    end 

end

%% Dump the paired comparisons into a table
all_paired = [];

for i = 1:length(paired.post_emu_def)
    all_paired = [all_paired;paired.post_emu_def(i).text];
    for j = 1:length(paired.post_emu_def(i).group)
        all_paired = [all_paired;paired.post_emu_def(i).group(j).text];
    end
end

%% Turn this into a table and save it as a csv file
outT2 = cell2table(all_paired,'VariableNames',...
    {'Variable','Pre-EMU','Post-EMU','Statistic','p-value'});
writetable(outT2,[results_folder,'paired_table.csv']);

%% Make a paired plot
figure
set(gcf,'position',[10 10 1001 350])
tiledlayout(1,2,'tilespacing','tight','padding','tight')
for i = 1:2
    nexttile
    numbers = paired.post_emu_def(i).group(1).numbers;
    pre_intervention = paired.post_emu_def(i).group(1).pre_intervention;
    post_intervention = paired.post_emu_def(i).group(1).post_intervention;

    % Paired plot
    pre_p = plot(numbers(pre_intervention,1),numbers(pre_intervention,2),'o','linewidth',2); % one color for pre-intervention
    hold on
    post_p = plot(numbers(post_intervention,1),numbers(post_intervention,2),'o','linewidth',2);

    max_all = max(numbers,[],'all');
    min_all = min(numbers,[],'all');
    plot([min_all max_all],[min_all max_all],'k--','linewidth',2)
    xlim([min_all max_all])
    ylim([min_all max_all])
    xlabel('Pre-EMU event frequency')
    ylabel('Post-EMU event frequency')
    set(gca,'fontsize',15)
    legend([pre_p,post_p],{'Pre-intervention','Post-intervention'},'fontsize',15)
    title(sprintf('%s',formatted_p_values(paired.post_emu_def(i).group(1).p)))
end
print(gcf,[results_folder,'paired_plots'],'-dpng')
