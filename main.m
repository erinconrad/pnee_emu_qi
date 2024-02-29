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
file_name = 'PNEEEMUPathway-Erincleaner_DATA_2024-02-29_0553.csv';
label_file_name = 'PNEEEMUPathway-Erincleaner_DATA_LABELS_2024-02-29_0553.csv';

% load redcap output file into table
T = readtable([data_folder,file_name],'ReadVariableNames',true);
lT = readtable([data_folder,label_file_name]);

% Get table variable names
var_names_T = T.Properties.VariableNames;
var_names_lT = lT.Properties.VariableNames;

%% Prep text results file
fname = [results_folder,'other_results.html'];
fid = fopen(fname,'a');

%% Identify pre-intervention and post-intervention patients
all_pts = ones(length(T.phase),1); all_pts = logical(all_pts);
pre = T.phase == 1;
post = T.phase == 2;

%% Define table variables
table1_vars = {'age','gender','known_to_have_pnee_prior',...
    'dual_diagnosis','follows_with_psych_pre','asm_pre','any_ed_pre','time_to_diagnosis_in_days'};
table2_vars = {'psychiatry_consult_obtaine',...
    'talked_to_patient',...
    'follow_up_with_psychiatris',...
    'did_patient_follow_up_with',...
    };
table3_vars = {'improvement_fifty',...
    'fiftyper_12mo_combined',...
    'qol',...
    'qol_v2',...
    'ed_visits_or_hospitalizati'};
outcome_predictors = {'psychiatry_consult_obtaine',...
    'follow_up_with_psychiatris',...
    'patient_understands_diagno',...
    'patient_agrees_with_diagno',...
   'asm_dc',...
   'time_to_diagnosis_2_years',...
   'fnd_psych'};
outcomes = {'improvement_fifty','fiftyper_12mo_combined','ed_visits_or_hospitalizati'};
supp_table_vars = {'work','home','social','private','relationships'};



%% Group variables according to type
% This will let me know what analysis to run
% If you edit the 2nd column, that will change the name of the final table
% variable
continuous_normalish = {'age','Age';...
    'since_discharge_event_freq','Event frequency since discharge';...
    'since_discharge_event_seve','Event severity since discharge';...
    'work','Work';...
    'home','Home management';...
    'social','Social activities';...
    'private','Private activities';...
    'relationships','Ability to maintain close relationships';...
    'qol','Quality of life (3 months)';...
    'qol_v2','Quality of life (12-20 months)'};
continuous_skewed = {'event_frequency_pre','Event frequency pre-EMU';...
    'follow_up_interval','Follow up interval';...
    'event_frequency_post','Event frequency 3 months post-EMU';...
    'time_to_diagnosis_in_days','Time to diagnosis in days';...
    'post_emu_freq_12mo_comb','Event frequency 12-20 months-post EMU'};
binary = {'known_to_have_pnee_prior','Previously known to have PNEE','all';...
    'dual_diagnosis','Dual diagnosis','all';...
    'follows_with_psych_pre', 'Followed with psychiatry before EMU','all';...
    'asm_pre','ASMs before EMU','all';...
    'asm_dc','ASMs discontinued or lowered on discharge','all';...
    'were_asm_for_an_indication','Non-epilepsy indication for ASMs','all';...
    'talked_to_patient','Was the patient reached by phone call','all';...
    'patient_arranged_psych','Did the patient arrange psych follow-up','all';...
    'study_team_member_arranged','Did the study team member arrange psych follow-up','all';...
    'did_patient_follow_up_with','Did the patient follow up with neurology','all';...
    'was_follow_up_scheduled','Was neurology follow-up scheduled','all';...
    'patient_understands_diagno','Does patient understand diagnosis of PNEE','followup';...
    'did_your_understanding_of','Did patient understanding of PNEE improve after EMU','followup';...
    'patient_agrees_with_diagno','Does patient agree with diagnosis of PNEE','followup';...
    'ed_visits_or_hospitalizati','Any ED visits or hospitalizations since discharge','followup';...
    'patient_scheduled_appointm','Did the patient schedule appointment for second opinion','followup';...
    'anti_seizure_medication_wa','Were ASMs changed or stopped in the EMU','followup';...
    'psychiatry_consult_obtaine','Was psychiatry consult obtained in the EMU','followup';...
    'follow_up_with_psychiatris','Did the patient schedule or complete psychiatry follow up','followup';...
    'psych_med_was_changes_star','Was a psychiatric medication changed or started','followup';...
    'improvement_fifty','>50% Improvement in Event Frequency? (3 months)','followup';...
    'any_ed_pre','Any ED visits year prior to admission?','all';...
    'time_to_diagnosis_2_years','Time to diagnosis >2 years?','followup';...
    'twelve_month_50_improvemen','>50% Improvement in Event frequency? (12 months)','all';...
    'fiftyper_12mo_combined','>50% Improvement in Event frequency? (12-20 months)','all';...
    'fnd_psych','Is outpatient therapy specifically addressing PNEE?','all'};
non_binary = {'gender','Gender'};


%% Prepare structure that will contain stats for every variable

% Start a running count of all variables
var_count = 0;

%% Fill up continuous normal comparisons in structure
% For now I'll do non parametric stats for consistency, but we can change
% this
ncontinuous_normal = size(continuous_normalish,1);

% Loop over the continuous normal variables
for i = 1:ncontinuous_normal
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = continuous_normalish{i,1};
    out.var(var_count).type = 'continuous normalish';
    out.var(var_count).pretty_name = continuous_normalish{i,2};

    % Get pre and post for this variable
    var = T.(continuous_normalish{i,1});
    var_pre = var(pre);
    var_post = var(post);

    % Get median and IQR for pre and post
    pre_median = nanmedian(var_pre);
    pre_iqr = prctile(var_pre,[25 75]);
    post_median = nanmedian(var_post);
    post_iqr = prctile(var_post,[25 75]);
    all_median = nanmedian(var);
    all_iqr = prctile(var,[25 75]);

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
    col_num = strcmp(var_names_T,continuous_normalish{i,1});

    % Get the corresponding variable name from the label table
    lt_variable = var_names_lT{col_num};

    % Prepare text
    stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',all_median,all_iqr(1),all_iqr(2),sum(~isnan(var))),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(var_pre))),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(var_post))),...
        sprintf('U: %1.1f',U),...
        sprintf('%s',formatted_p_values(p))};

    % Add this info to the struct
    out.var(var_count).all_median = all_median;
    out.var(var_count).all_iqr = all_iqr;
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
ncontinuous_skewed = size(continuous_skewed,1);

% Loop over the continuous normal variables
for i = 1:ncontinuous_skewed
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = continuous_skewed{i,1};
    out.var(var_count).type = 'continuous skewed';
    out.var(var_count).pretty_name = continuous_skewed{i,2};

    % Get pre and post for this variable
    var = T.(continuous_skewed{i,1});
    var_pre = var(pre);
    var_post = var(post);

    % Get median and IQR for pre and post
    all_median = nanmedian(var);
    all_iqr = prctile(var,[25 75]);
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
    col_num = strcmp(var_names_T,continuous_skewed{i,1});

    % Get the corresponding variable name from the label table
    lt_variable = var_names_lT{col_num};

    %{
    % Prepare text
    stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
        sprintf('%1.1f (%1.1f-%1.1f)',all_median,all_iqr(1),all_iqr(2)),...
        sprintf('%1.1f (%1.1f-%1.1f)',pre_median,pre_iqr(1),pre_iqr(2)),...
        sprintf('%1.1f (%1.1f-%1.1f)',post_median,post_iqr(1),post_iqr(2)),...
        sprintf('U: %1.1f',U),...
        sprintf('%s',formatted_p_values(p))};
    %}
    % Prepare text
    stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',all_median,all_iqr(1),all_iqr(2),sum(~isnan(var))),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(var_pre))),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(var_post))),...
        sprintf('U: %1.1f',U),...
        sprintf('%s',formatted_p_values(p))};


    % Add this info to the struct
    out.var(var_count).all_median = all_median;
    out.var(var_count).all_iqr = all_iqr;
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
nbinary = size(binary,1);

% Loop over the binary variables
for i = 1:nbinary
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = binary{i,1};
    out.var(var_count).type = 'binary';
    out.var(var_count).pretty_name = binary{i,2};

    

    % Get pre and post for this variable
    var = T.(binary{i,1});
    var_pre = var(pre);
    var_post = var(post);

    % fix for talked to patient
    if strcmp(binary{i,1},'talked_to_patient')
        var_pre(isnan(var_pre)) = 0;
    end

    % Get the column number of this table
    col_num = strcmp(var_names_T,binary{i,1});

    % Get the corresponding column from the label table
    lt_variable = var_names_lT{col_num};
    lT_col = lT.(lt_variable);

    % double check that it really is binary
    assert(sum(var == 1 | var == 0 | isnan(var)) == length(var))

    % get the labels corresponding to 1 and 0
    label1 = unique(lT_col(var==1));
    label0 = unique(lT_col(var==0));
    assert(length(label1)<=1); assert(length(label0)<=1);

    if ~iscell(label1)
        if label1 == 1 && label0 == 0
            label1 = {'Yes'}; label0 = {'No'};
        end
    end

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
    all1 = sum(var==1);
    all0 = sum(var==0);

    % Prep 2x2 table
    tbl_2x2 = table([pre0;pre1],[post0;post1],...
        'VariableNames',{'Pre','Post'},'RowNames',{'0','1'});
    [~,p,stats] = fishertest(tbl_2x2);

    % Prep text
    %{
    stats_text = {sprintf('%s',out.var(var_count).pretty_name),'','','',...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label0{1}),...
        sprintf('%d (%1.1f%%)',all0,all0/length(var)*100),...
        sprintf('%d (%1.1f%%)',pre0,pre0/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post0,post0/length(var_post)*100),'','';...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',all1,all1/length(var)*100),...
        sprintf('%d (%1.1f%%)',pre1,pre1/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post1,post1/length(var_post)*100),'',''};
    %}
    stats_text = {sprintf('%s',out.var(var_count).pretty_name),...
        sprintf('N = %d',sum(~isnan(var))),...
        sprintf('N = %d',sum(~isnan(var_pre))),...
        sprintf('N = %d',sum(~isnan(var_post))),...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label0{1}),...
        sprintf('%d (%1.1f%%)',all0,all0/(all0+all1)*100),...
        sprintf('%d (%1.1f%%)',pre0,pre0/(pre0+pre1)*100),...
        sprintf('%d (%1.1f%%)',post0,post0/(post0+post1)*100),'','';...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',all1,all1/(all0+all1)*100),...
        sprintf('%d (%1.1f%%)',pre1,pre1/(pre0+pre1)*100),...
        sprintf('%d (%1.1f%%)',post1,post1/(post0+post1)*100),'',''};

    % Add this to structure
    out.var(var_count).tbl = tbl_2x2;
    out.var(var_count).p = p;
    out.var(var_count).stats = stats;
    out.var(var_count).text = stats_text;


end

%% Fill up non-binary comparisons
nnon_binary = size(non_binary,1);
% Loop over the non binary variables
for i = 1:nnon_binary
    var_count = var_count + 1; % increase the index of the variable

    % Get variable name
    out.var(var_count).name = non_binary{i,1};
    out.var(var_count).type = 'non-binary';
    out.var(var_count).pretty_name = non_binary{i,2};

    % Get pre and post for this variable
    var = T.(non_binary{i,1});
    var_pre = var(pre);
    var_post = var(post);

    % double check that options are just 1, 2, and 3
    assert(sum(var <= 1 | var <= 2 | var <= 3) == length(var))

    % Get the column number of this table
    col_num = strcmp(var_names_T,non_binary{i,1});

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
    all1 = sum(var==1);
    all2 = sum(var==2);
    all3 = sum(var==3);

    % Prep 2x3 table
    tbl_2x3 = table([pre1;pre2;pre3],[post1;post2;post3],...
        'VariableNames',{'Pre','Post'},'RowNames',{'1','2','3'});
    
    % Convert this into the format needed for MyFisher23
    x = table2array(tbl_2x3);
    x = x';

    % do myfisher23
    p = myfisher23(x);

    % Prep text
    stats_text = {sprintf('%s',out.var(var_count).pretty_name),sprintf('N = %d',sum(~isnan(var))),...
        sprintf('N = %d',sum(~isnan(var_pre))),...
        sprintf('N = %d',sum(~isnan(var_post))),'',...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',all1,all1/length(var)*100),...
        sprintf('%d (%1.1f%%)',pre1,pre1/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post1,post1/length(var_post)*100),'','';...
        sprintf('%s',label2{1}),...
        sprintf('%d (%1.1f%%)',all2,all2/length(var)*100),...
        sprintf('%d (%1.1f%%)',pre2,pre2/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post2,post2/length(var_post)*100),'','';...
        sprintf('%s',label3{1}),...
        sprintf('%d (%1.1f%%)',all3,all3/length(var)*100),...
        sprintf('%d (%1.1f%%)',pre3,pre3/length(var_pre)*100),...
        sprintf('%d (%1.1f%%)',post3,post3/length(var_post)*100),'',''};

    % Add this to structure
    out.var(var_count).tbl = tbl_2x3;
    out.var(var_count).p = p;
    out.var(var_count).text = stats_text;

end




%% Look at pre-post for patients who saw or scheduled psych
%{
var_count = var_count + 1;
saw_or_scheduled_psych = T.follow_up_with_psychiatris;
out.var(var_count).name = 'improvement_pp';
out.var(var_count).type = 'continuous normalish';
out.var(var_count).pretty_name = 'Improvement in frequency for those who scheduled or saw psych';
var_pre = T.since_discharge_event_freq(pre);
var_post = T.since_discharge_event_freq(post & saw_or_scheduled_psych==1); % re-define post to also require seeing or scheduling psych

% Get median and IQR for pre and post
pre_median = nanmedian(var_pre);
pre_iqr = prctile(var_pre,[25 75]);
post_median = nanmedian(var_post);
post_iqr = prctile(var_post,[25 75]);

[p,~,stats] = ranksum(var_pre,var_post);
W = stats.ranksum;

% Get U statistic. DOUBLE CHECK THAT THIS IS RIGHT
U = W-sum(~isnan(var_pre))*(sum(~isnan(var_pre))+1)/2;

% Prepare text
stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
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
%}


%% Dump all the pre vs post comparisons into one huge table
all_pre_vs_post_comp = [];
n_pre_vs_post_comp = length(out.var);
for i = 1:n_pre_vs_post_comp
    all_pre_vs_post_comp = [all_pre_vs_post_comp;...
        out.var(i).text];
end

%% Turn this into a table and save it as a csv file
outT = cell2table(all_pre_vs_post_comp,'VariableNames',...
    {'Variable','All','Pre','Post','Statistic','p-value'});
writetable(outT,[results_folder,'pre_vs_post_table.csv']);

%% Make Table 1
% find table 1 vars
all_table1_comp = [];
for i = 1:length(table1_vars)
    for j = 1:length(out.var)
        if strcmp(out.var(j).name,table1_vars{i})
            all_table1_comp = [all_table1_comp;out.var(j).text];
        end
    end
end
table1T = cell2table(all_table1_comp,'VariableNames',...
    {'Variable','All','Pre','Post','Statistic','p-value'});
writetable(table1T,[results_folder,'table1.csv'],'QuoteStrings',true);

%% Make Table 2
% find table 2 vars
all_table2_comp = [];
for i = 1:length(table2_vars)
    for j = 1:length(out.var)
        if strcmp(out.var(j).name,table2_vars{i})
            all_table2_comp = [all_table2_comp;out.var(j).text];
        end
    end
end
table2T = cell2table(all_table2_comp,'VariableNames',...
    {'Variable','All','Pre','Post','Statistic','p-value'});
writetable(table2T,[results_folder,'table2.csv'],'QuoteStrings',true);

%% Make Table 3
% find table 3 vars
all_table3_comp = [];
for i = 1:length(table3_vars)
    for j = 1:length(out.var)
        if strcmp(out.var(j).name,table3_vars{i})
            all_table3_comp = [all_table3_comp;out.var(j).text];
        end
    end
end
table3T = cell2table(all_table3_comp,'VariableNames',...
    {'Variable','All','Pre','Post','Statistic','p-value'});
writetable(table3T,[results_folder,'table3.csv'],'QuoteStrings',true);

%% Table 4 - outcome predictors
outcome_pred_text = cell(length(outcome_predictors),1);
for i = 1:length(outcome_predictors)
    for j = 1:length(outcomes)
        curr = T.(outcome_predictors{i});
        curr_out = T.(outcomes{j});
    
        nonan = ~any([isnan(curr),isnan(curr_out)],2);
        curr_nonan = curr(nonan);
        curr_out_nonan = curr_out(nonan);
    
        if 0
            table(curr,curr_out)
        end
    
        tbl_2x2 = crosstab(curr_nonan,curr_out_nonan);
        [~,p,stats] = fishertest(tbl_2x2);
    
        for k = 1:length(out.var)
            if strcmp(out.var(k).name,outcome_predictors{i})
                curr_name = out.var(k).pretty_name;
            end
        end
        if j == 1
            outcome_pred_text{i} = [outcome_pred_text{i},...
                {curr_name,sprintf('%1.2f (%1.2f-%1.2f)\nN=%d',stats.OddsRatio,...
                stats.ConfidenceInterval(1),stats.ConfidenceInterval(2),sum(nonan)),...
                sprintf('%1.2f',p)}];
        else
            outcome_pred_text{i} = [outcome_pred_text{i},...
                {sprintf('%1.2f (%1.2f-%1.2f)\nN=%d',stats.OddsRatio,...
                stats.ConfidenceInterval(1),stats.ConfidenceInterval(2),sum(nonan)),...
                sprintf('%1.2f',p)}];
        end
    end
end

combined_array = outcome_pred_text{1};
for i = 2:length(outcome_pred_text)
    combined_array = vertcat(combined_array, outcome_pred_text{i});
end

table4T = cell2table(combined_array,'VariableNames',...
    {'Variable','OR (95% CI) for >50% event reduction (3 months)','p-value1',...
    'OR (95% CI) for >50% event reduction (12-20 months)','p-value2',...
    'ED visits or hospitalizations','p-value3'});
writetable(table4T,[results_folder,'table4.csv'],'QuoteStrings',true);

%% Table 5
table5_text = cell(2,1);
pre_freq_var = 'event_frequency_pre';
post_freq_vars = {'event_frequency_post','post_emu_freq_12mo_comb'};
npost = length(post_freq_vars);

% Loop over two ways of defining post-EMU event frequency (visit vs survey)
for iphase = 1:2
    if iphase == 1
        group = pre;
        which_group = 'Care-as-usual';
    else
        group = post;
        which_group = 'Intervention';
    end

    for i = 1:npost
    
        % Get pre- and post-EMU 
        pre_emu = T.(pre_freq_var);
        post_emu = T.(post_freq_vars{i});
    
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
        if i == 1
            table5_text{iphase} = [table5_text{iphase},{sprintf('%s',which_group),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(pre_emu(group)))),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(post_emu(group)))),...
            sprintf('T+: %1.1f',signedrank),...
            sprintf('%s',formatted_p_values(p))}];
        else
            table5_text{iphase} = [table5_text{iphase},{sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(post_emu(group)))),...
            sprintf('T+: %1.1f',signedrank),...
            sprintf('%s',formatted_p_values(p))}];
        end

    end 

end

combined_array = table5_text{1};
for i = 2:length(table5_text)
    combined_array = vertcat(combined_array, table5_text{i});
end

table5T = cell2table(combined_array,'VariableNames',...
    {'Phase','Pre-EMU events/month','3 months post-EMU events/month',...
    'Statistic1','p-value1',...
    '12-20 months post-EMU events/month','Statistic2','p-value2'});
writetable(table5T,[results_folder,'table5.csv'],'QuoteStrings',true);


%% Supp table 1
all_supp_table_comp = [];
for i = 1:length(supp_table_vars)
    for j = 1:length(out.var)
        if strcmp(out.var(j).name,supp_table_vars{i})
            all_supp_table_comp = [all_supp_table_comp;out.var(j).text];
        end
    end
end
tableST = cell2table(all_supp_table_comp,'VariableNames',...
    {'Variable','All','Pre','Post','Statistic','p-value'});
writetable(tableST,[results_folder,'supp_table.csv'],'QuoteStrings',true);

%% Prepare paired analyses
pre_freq_var = 'event_frequency_pre';
post_freq_vars = {'event_frequency_post','post_emu_freq_12mo_comb'};
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
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(pre_emu(group)))),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(post_emu(group)))),...
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



if 0
    figure
    nexttile
    histogram(T.twelve_month_relative_chan)
    nexttile
    histogram(T.twelve_month_phone_relative_chan)
end

%% Turn this into a table and save it as a csv file
outT2 = cell2table(all_paired,'VariableNames',...
    {'Variable','Pre-EMU','Post-EMU','Statistic','p-value'});
writetable(outT2,[results_folder,'paired_table.csv'],'QuoteStrings',true);

%% Make a paired plot
if 1
figure
set(gcf,'position',[10 10 1400 1000])
tiledlayout(2,2,'tilespacing','compact','padding','tight')
for ip = 1:2
numbers = paired.post_emu_def(ip).group(1).numbers;
for i = 1:2
    nexttile
    if i == 1
        
        arm_text = 'Care-as-usual';
        data = numbers(paired.post_emu_def(ip).group(1).pre_intervention,:);
        p =  paired.post_emu_def(ip).group(2).p;
    else
        arm_text = 'Intervention';
        data = numbers(paired.post_emu_def(ip).group(1).post_intervention,:);
        p =  paired.post_emu_def(ip).group(3).p;
    end

    % paired plot
    plot([1 2],data,'-k','linewidth',2)
    hold on
    xlim([0.5 2.5])
    xticks([1 2])
    if ip == 1
        xticklabels({'Pre-EMU','3 months post-EMU'})
    else
        xticklabels({'Pre-EMU','12-20 months post-EMU'})
    end
    title(arm_text)
    if i == 1
        ylabel('Events/month')
    end
    set(gca,'fontsize',25)

    yl = ylim;
    ybar_top = yl(1)+(yl(2)-yl(1))*1.08;
    ybar_ticks = yl(1)+(yl(2)-yl(1))*1.05;
    ytext = yl(1)+(yl(2)-yl(1))*1.15;
    yl_new = [yl(1) yl(1)+(yl(2)-yl(1))*1.25];
    plot([1 2],[ybar_top ybar_top],'k-','linewidth',2)
    plot([1 1],[ybar_ticks ybar_top],'k-','linewidth',2)
    plot([2 2],[ybar_ticks ybar_top],'k-','linewidth',2)
    if p < 0.001
        ptext = 'p < 0.001';
    elseif p < 0.05
        ptext = sprintf('p = %1.3f',p);
    else
        ptext = sprintf('p = %1.2f',p);
    end
    text(1.5,ytext,ptext,'fontsize',25,'HorizontalAlignment','center')
    ylim(yl_new)
    

end
end
set(gcf,'renderer','painters')
print(gcf,[results_folder,'paired_plots2'],'-dpng')
end


if 0
figure
set(gcf,'position',[10 10 550 400])
tiledlayout(1,1,'tilespacing','tight','padding','tight')
for i = 1
    nexttile
    numbers = paired.post_emu_def(i).group(1).numbers;
    pre_intervention = paired.post_emu_def(i).group(1).pre_intervention;
    post_intervention = paired.post_emu_def(i).group(1).post_intervention;

    % Paired plot
    pre_p = plot(numbers(pre_intervention,1),numbers(pre_intervention,2),'o','linewidth',2,...
        'markersize',17); % one color for pre-intervention
    hold on
    post_p = plot(numbers(post_intervention,1),numbers(post_intervention,2),'x','linewidth',2,...
        'markersize',17);

    max_all = max(numbers,[],'all');
    min_all = min(numbers,[],'all');
    plot([min_all max_all],[min_all max_all],'k--','linewidth',2)
    xlim([min_all max_all])
    ylim([min_all max_all])
    xlabel('Pre-EMU events/month')
    ylabel('Post-EMU events/month')
    set(gca,'fontsize',15)
    legend([pre_p,post_p],{'Pre-intervention','Post-intervention'},'fontsize',25,...
        'location','northwest')
    set(gca,'fontsize',25)
    %title(sprintf('%s',formatted_p_values(paired.post_emu_def(i).group(1).p)))
end
set(gcf,'renderer','painters')

print(gcf,[results_folder,'paired_plots'],'-dpng')
end

%% Additional analysis
%{
% Percentage of patients without prior psych care who arrange appt
no_prior_psych_arrange_appt_all = T.follows_with_psych_pre == 0 & (T.patient_arranged_psych==1 | T.study_team_member_arranged == 1);
no_prior_psych_arrange_appt_pt = T.follows_with_psych_pre == 0 & T.patient_arranged_psych==1;
no_prior_psych_arrange_appt_team = T.follows_with_psych_pre == 0 & T.study_team_member_arranged==1;
fprintf(fid,'<p>Number of patients without prior psych care: %d.</p>',sum(T.follows_with_psych_pre == 0));
fprintf(fid,['<p>Number (%%) of these who arranged psych follow up: %d (%1.1f%%). '...
    'Of these, %d (%1.1f%%) arranged follow up themselves, and %d (%1.1f%%) had follow-up '...
    'arranged by team.</p>'],sum(no_prior_psych_arrange_appt_pt),sum(no_prior_psych_arrange_appt_pt)/sum(T.follows_with_psych_pre == 0)*100,...
    sum(no_prior_psych_arrange_appt_pt),sum(no_prior_psych_arrange_appt_pt)/sum(no_prior_psych_arrange_appt_all)*100,...
    sum(no_prior_psych_arrange_appt_team),sum(no_prior_psych_arrange_appt_team)/sum(no_prior_psych_arrange_appt_all)*100);

% Get no psych follow up reasons
no_psych_fields = var_names_lT(contains(var_names_lT,'IfNoPsychFollo'));
assert(length(no_psych_fields)==20)% confirm it's 20
no_psych_fields = no_psych_fields(1:10); % take the first 10
no_psych_text = cell(10,1);

no_psych_numbers = nan(10,1);
no_psych_numbers_pre = nan(10,1);
no_psych_numbers_post = nan(10,1);
for i = 1:10
    curr = no_psych_fields{i};
    C = strsplit(curr,'_'); 
    no_psych_text{i} = C{end};
    no_psych_numbers(i) = sum(strcmp(lT.(curr),'Checked'));
    no_psych_numbers_pre(i) = sum(strcmp(lT.(curr),'Checked') & strcmp(lT.Phase,'Pre-Intervention'));
    no_psych_numbers_post(i) = sum(strcmp(lT.(curr),'Checked') & strcmp(lT.Phase,'Post-Intervention'));
end
no_psych_text{3} = 'Insurance';
no_psych_text{9} = 'Scheduled';
comb = [no_psych_numbers_pre';no_psych_numbers_post'];
if 0
figure
bar(1:10,comb)
xticks(1:size(comb,2))
xticklabels(no_psych_text)
legend({'Pre','Post'})
end
%}

%{
Old fields
continuous_normalish = {'age','Age';...
    'since_discharge_event_freq','Event frequency since discharge';...
    'since_discharge_event_seve','Event severity since discharge';...
    'work','Work';...
    'home','Home';...
    'social','Social';...
    'private','Private';...
    'relationships','Relationships';...
    'qol','Quality of life';...
    'work1','Work (survey)';...
    'home1','Home (survey)';...
    'social1','Social (survey)';...
    'private1','Private (survey)';...
    'relationships1','Relationships (survey)';...
    'qol1','Quality of life (survey)';...
    'since_discharge_event_freq1','Event frequency since discharge (survey)';...
    'since_discharge_event_seve1','Event severity since discharge (survey)'};
continuous_skewed = {'event_frequency_pre','Event frequency pre-EMU';...
    'follow_up_interval','Follow up interval';...
    'event_frequency_post','Event frequency post-EMU';...
    'event_frequency_post1','Event frequency post-EMU (survey)';...
    'time_to_diagnosis_in_days','Time to diagnosis in days';...
    'event_frequency_post1_2b1583','Event frequency 12 months-post EMU'};
binary = {'known_to_have_pnee_prior','Previously known to have PNEE','all';...
    'dual_diagnosis','Dual diagnosis','all';...
    'follows_with_psych_pre', 'Followed with psychiatry before EMU','all';...
    'asm_pre','ASMs before EMU','all';...
    'asm_dc','ASMs discontinued or lowered on discharge','all';...
    'were_asm_for_an_indication','Non-epilepsy indication for ASMs','all';...
    'talked_to_patient','Was the patient reached by phone call','all';...
    'patient_arranged_psych','Did the patient arrange psych follow-up','all';...
    'study_team_member_arranged','Did the study team member arrange psych follow-up','all';...
    'did_patient_follow_up_with','Did the patient follow up with neurology','all';...
    'was_follow_up_scheduled','Was neurology follow-up scheduled','all';...
    'patient_understands_diagno','Does patient understand diagnosis of PNEE','followup';...
    'did_your_understanding_of','Did patient understanding of PNEE improve after EMU','followup';...
    'patient_agrees_with_diagno','Does patient agree with diagnosis of PNEE','followup';...
    'ed_visits_or_hospitalizati','Any ED visits or hospitalizations since discharge','followup';...
    'patient_scheduled_appointm','Did the patient schedule appointment for second opinion','followup';...
    'anti_seizure_medication_wa','Were ASMs changed or stopped in the EMU','followup';...
    'psychiatry_consult_obtaine','Was psychiatry consult obtained in the EMU','followup';...
    'follow_up_with_psychiatris','Did the patient schedule or complete psychiatry follow up','followup';...
    'psych_med_was_changes_star','Was a psychiatric medication changed or started','followup';...
    'self_help_or_apps_used_for','Were self-help or apps used for PNEE','followup';...
    'patient_understands_diagno1','Does patient understand diagnosis of PNEE (survey)','all';...
    'did_your_understanding_of1','Did patient understanding of PNEE improve after EMU (survey)','all';...
    'patient_agrees_with_diagno1','Does patient agree with diagnosis of PNEE (survey)','all';...
    'ed_visits_or_hospitalizati1','Any ED visits or hospitalizations since discharge (survey)','all';...
    'patient_scheduled_appointm1','Did the patient schedule appointment for second opinion (survey)','all';...
    'anti_seizure_medication_wa1','Were ASMs changed or stopped in the EMU (survey)','all';...
    'follow_up_with_psychiatris1','Did the patient schedule or complete psychiatry follow up (survey)','all';...
    'psych_med_was_changes_star1','Was a psychiatric medication changed or started (survey)','all';...
    'self_help_or_apps_used_for1','Were self-help or apps used for PNEE (survey)','all';...
    'improvement_fifty','>50% Improvement in Event Frequency? (3 months)','followup';...
    'any_ed_pre','Any ED visits year prior to admission?','all';...
    'time_to_diagnosis_2_years','Time to diagnosis >2 years?','followup';...
    'twelve_month_50_improvemen','>50% Improvement in Event frequency? (12 months)','all'};
non_binary = {'gender','Gender'};

%}