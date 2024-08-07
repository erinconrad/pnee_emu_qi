clear
clc


%% Dependencies
% This requires adding myFisher23 to your path 
% https://www.mathworks.com/matlabcentral/fileexchange/15399-myfisher23

%% Get file locations and load redcap file
locations = pnee_locations;
data_folder = locations.data;
results_folder = locations.results;
file_name = 'pnee_data.csv';
label_file_name = 'pnee_labels.csv';

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
table1_vars = {'age','gender','event_frequency_pre','known_to_have_pnee_prior',...
    'dual_diagnosis','follows_with_psych_pre','asm_pre','any_ed_pre','time_to_diagnosis_in_days'};
table2_vars = {'psychiatry_consult_obtaine',...
    'talked_to_patient',...
    'follow_up_with_psychiatris'};
table3_vars = {'rel_diff_pre_post',...
    'resolution_at_3_months',...
    'rel_change_12_comb',...
    'resolution_combined_12',...
    'follow_up_with_psychiatris_v2',...
    'comb_ed_or_hospital_visits',...
    'qol_v2'};
outcome_predictors = {'psychiatry_consult_obtaine',...
    'follow_up_with_psychiatris',...
    'follow_up_with_psychiatris_v2',...
    'fnd_psych',...
    'patient_understands_diagno',...
    'patient_agrees_with_diagno',...
   'asm_dc',...
   'time_to_diagnosis_2_years'};
outcomes = {'resolution_combined_12','comb_ed_or_hospital_visits'};
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
    'post_emu_freq_12mo_comb','Event frequency 12-20 months-post EMU';...
    'rel_diff_pre_post','Percent reduction in event frequency (3 months)';...
    'rel_change_12_comb','Percent reduction in event frequency (12-20 months)'};
binary = {'known_to_have_pnee_prior','Previously known to have FS/PNEE','all';...
    'dual_diagnosis','Dual diagnosis','all';...
    'follows_with_psych_pre', 'Followed with psych before EMU','all';...
    'asm_pre','ASMs before EMU','all';...
    'asm_dc','ASMs discontinued or lowered on discharge','all';...
    'talked_to_patient','Was the patient reached by phone call','all';...
    'patient_arranged_psych','Did the patient arrange psych follow-up','all';...
    'study_team_member_arranged','Did the study team member arrange psych follow-up','all';...
    'did_patient_follow_up_with','Did the patient follow up with neurology','all';...
    'was_follow_up_scheduled','Was neurology follow-up scheduled','all';...
    'patient_understands_diagno','Does patient understand diagnosis of PNEE','followup';...
    'did_your_understanding_of','Did patient understanding of PNEE improve after EMU','followup';...
    'patient_agrees_with_diagno','Does patient agree with diagnosis of PNEE','followup';...
    'comb_ed_or_hospital_visits','Any known ED visits or hospitalizations since discharge','all';...
    'patient_scheduled_appointm','Did the patient schedule appointment for second opinion','followup';...
    'anti_seizure_medication_wa','Were ASMs changed or stopped in the EMU','followup';...
    'psychiatry_consult_obtaine','Was psychiatry consult obtained in the EMU','followup';...
    'follow_up_with_psychiatris','Did the patient schedule psychiatry follow up within 3 months?','followup';...
    'psych_med_was_changes_star','Was a psychiatric medication changed or started','followup';...
    'improvement_fifty','>50% Improvement in Event Frequency? (3 months)','followup';...
    'any_ed_pre','Any ED visits year prior to admission?','all';...
    'time_to_diagnosis_2_years','Time to diagnosis >2 years?','followup';...
    'twelve_month_50_improvemen','>50% Improvement in Event frequency? (12 months)','all';...
    'fiftyper_12mo_combined','>50% Improvement in Event frequency? (12-20 months)','all';...
    'fnd_psych','Is outpatient therapy specifically addressing PNEE?','all';...
    'resolution_combined_12','Event resolution (12-20 months)','all';...
    'resolution_at_3_months','Event resolution (3 months)','all';...
    'follow_up_with_psychiatris_v2','Did the patient complete psychiatry follow up within 12-20 months?','all'};
non_binary = {'gender','Gender'};
% 'were_asm_for_an_indication','Non-epilepsy indication for ASMs','all';...

%% Prepare structure that will contain stats for every variable

% Start a running count of all variables
var_count = 0;

%% Say people who have no info for whether meds reduced did NOT have them reduced if they were not on meds to begin with
no_meds_initially = strcmp(lT.OnSeizureMedsPriorToAdmission_,'No');
T.asm_dc(no_meds_initially) = 0;
lT.IfYes_WereTheyDiscontinuedOrLowered_(no_meds_initially) = {'No'};

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

    if strcmp(out.var(var_count).name,'rel_diff_pre_post') || strcmp(out.var(var_count).name,'rel_change_12_comb')

        % Prepare text
        stats_text_missing = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)\n(N = %d of %d with data)',all_median,all_iqr(1),all_iqr(2),sum(~isnan(var)),length(var)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)\n(N = %d of %d with data)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(var_pre)),length(var_pre)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)\n(N = %d of %d with data)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(var_post)),length(var_post)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};
    
        stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)',all_median,all_iqr(1),all_iqr(2)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)',pre_median,pre_iqr(1),pre_iqr(2)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)',post_median,post_iqr(1),post_iqr(2)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};

    else

        % Prepare text
        stats_text_missing = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',all_median,all_iqr(1),all_iqr(2),sum(~isnan(var)),length(var)),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(var_pre)),length(var_pre)),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(var_post)),length(var_post)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};
    
        stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f (%1.1f-%1.1f)',all_median,all_iqr(1),all_iqr(2)),...
            sprintf('%1.1f (%1.1f-%1.1f)',pre_median,pre_iqr(1),pre_iqr(2)),...
            sprintf('%1.1f (%1.1f-%1.1f)',post_median,post_iqr(1),post_iqr(2)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};
    end

    if strcmp(out.var(var_count).name,'comb_ed_or_hospital_visits')
        out.var(var_count).text_missing = stats_text;
    else
        out.var(var_count).text_missing = stats_text_missing;
    end

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
    out.var(var_count).text_missing = stats_text_missing;

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
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',all_median,all_iqr(1),all_iqr(2),sum(~isnan(var))),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(var_pre))),...
        sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(var_post))),...
        sprintf('U: %1.1f',U),...
        sprintf('%s',formatted_p_values(p))};
    %}
    if strcmp(out.var(var_count).name,'rel_diff_pre_post') || strcmp(out.var(var_count).name,'rel_change_12_comb')

        % Prepare text
        stats_text_missing = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)\n(N = %d of %d with data)',all_median,all_iqr(1),all_iqr(2),sum(~isnan(var)),length(var)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)\n(N = %d of %d with data)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(var_pre)),length(var_pre)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)\n(N = %d of %d with data)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(var_post)),length(var_post)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};
    
        stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)',all_median,all_iqr(1),all_iqr(2)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)',pre_median,pre_iqr(1),pre_iqr(2)),...
            sprintf('%1.1f%% (%1.1f-%1.1f%%)',post_median,post_iqr(1),post_iqr(2)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};

    else

        stats_text_missing = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',all_median,all_iqr(1),all_iqr(2),sum(~isnan(var)),length(var)),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(var_pre)),length(var_pre)),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(var_post)),length(var_post)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};
    
        stats_text = {sprintf('%s: median (IQR)',out.var(var_count).pretty_name),...
            sprintf('%1.1f (%1.1f-%1.1f)',all_median,all_iqr(1),all_iqr(2)),...
            sprintf('%1.1f (%1.1f-%1.1f)',pre_median,pre_iqr(1),pre_iqr(2)),...
            sprintf('%1.1f (%1.1f-%1.1f)',post_median,post_iqr(1),post_iqr(2)),...
            sprintf('U: %1.1f',U),...
            sprintf('%s',formatted_p_values(p))};
    end


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
    out.var(var_count).text_missing = stats_text_missing;


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

    % Get the column number of this table
    col_num = strcmp(var_names_T,binary{i,1});

    % Get the corresponding column from the label table
    lt_variable = var_names_lT{col_num};
    lT_col = lT.(lt_variable);


    % fix for talked to patient
    if strcmp(binary{i,1},'talked_to_patient')
        var_pre(isnan(var_pre)) = 0;
        var(isnan(var)) = 0;
        lT_col(cellfun(@isempty,lT_col)) = {'No'};
    end

    
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
    prenan = sum(isnan(var_pre));
    postnan = sum(isnan(var_post));
    allnan = sum(isnan(var));

    assert(pre1+pre0+prenan == length(var_pre))
    assert(post1+post0+postnan == length(var_post))
    assert(all1+all0+allnan == length(var))

    % Prep 2x2 table
    tbl_2x2 = table([pre0;pre1],[post0;post1],...
        'VariableNames',{'Pre','Post'},'RowNames',{'0','1'});
    [~,p,stats] = fishertest(tbl_2x2);

    % Prep text
    %{
    stats_text = {sprintf('%s',out.var(var_count).pretty_name),...
        sprintf('N = %d',sum(~isnan(var))),...
        sprintf('N = %d',sum(~isnan(var_pre))),...
        sprintf('N = %d',sum(~isnan(var_post))),...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label0{1}),...
        sprintf('%d (%1.1f%%)',all0,all0/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre0,pre0/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post0,post0/(length(var_post))*100),'','';...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',all1,all1/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre1,pre1/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post1,post1/(length(var_post))*100),'',''};

    stats_text_missing = {sprintf('%s',out.var(var_count).pretty_name),...
        sprintf('N = %d',length(var)),...
        sprintf('N = %d',length(var_pre)),...
        sprintf('N = %d',length(var_post)),...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label0{1}),...
        sprintf('%d (%1.1f%%)',all0,all0/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre0,pre0/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post0,post0/(length(var_post))*100),'','';...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',all1,all1/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre1,pre1/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post1,post1/(length(var_post))*100),'','';...
        sprintf('Unknown'),...
        sprintf('%d (%1.1f%%)',allnan,allnan/(length(var))*100),...
        sprintf('%d (%1.1f%%)',prenan,prenan/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',postnan,postnan/(length(var_post))*100),'',''};
    %}

    stats_text = {sprintf('%s',out.var(var_count).pretty_name),...
        '',...
        '',...
        '',...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label0{1}),...
        sprintf('%d (%1.1f%%)',all0,all0/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre0,pre0/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post0,post0/(length(var_post))*100),'','';...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',all1,all1/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre1,pre1/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post1,post1/(length(var_post))*100),'',''};

    stats_text_missing = {sprintf('%s',out.var(var_count).pretty_name),...
        '',...
        '',...
        '',...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('%s',label0{1}),...
        sprintf('%d (%1.1f%%)',all0,all0/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre0,pre0/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post0,post0/(length(var_post))*100),'','';...
        sprintf('%s',label1{1}),...
        sprintf('%d (%1.1f%%)',all1,all1/(length(var))*100),...
        sprintf('%d (%1.1f%%)',pre1,pre1/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',post1,post1/(length(var_post))*100),'','';...
        sprintf('Unknown'),...
        sprintf('%d (%1.1f%%)',allnan,allnan/(length(var))*100),...
        sprintf('%d (%1.1f%%)',prenan,prenan/(length(var_pre))*100),...
        sprintf('%d (%1.1f%%)',postnan,postnan/(length(var_post))*100),'',''};

    % Add this to structure
    out.var(var_count).tbl = tbl_2x2;
    out.var(var_count).p = p;
    out.var(var_count).stats = stats;
    out.var(var_count).text = stats_text;
    out.var(var_count).text_missing = stats_text_missing;


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
    %{
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
    %}
    stats_text = {sprintf('%s',out.var(var_count).pretty_name),'',...
        '',...
        '','',...
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

%% Fig 1 info for Kelly
% 12-20 care as usual follow up via phone call
n_12_cau_ph = sum(~isnan(T.fiftyper_12mon_phonecall)&T.phase==1);
% 12-20 care as usual follow up via chart review
n_12_cau_ch = sum((isnan(T.fiftyper_12mon_phonecall) & ~isnan(T.event_frequency_post1_2b1583))&T.phase==1);
% 12-20 cau without phone call
n_12_cau_no_ph = sum(isnan(T.fiftyper_12mon_phonecall)&T.phase==1);

% 12-20 intervention follow up via phone call
n_12_int_ph = sum(~isnan(T.fiftyper_12mon_phonecall)&T.phase==2);
% 12-20 intervention follow up via chart review
n_12_int_ch = sum((isnan(T.fiftyper_12mon_phonecall) & ~isnan(T.event_frequency_post1_2b1583))&T.phase==2);
% 12-20 intervention no phone
n_12_int_no_ph = sum(isnan(T.fiftyper_12mon_phonecall)&T.phase==2);

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
            all_table2_comp = [all_table2_comp;out.var(j).text_missing];
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
            all_table3_comp = [all_table3_comp;out.var(j).text_missing];
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

        match = 0;
        for k = 1:length(out.var)
            if strcmp(out.var(k).name,outcome_predictors{i})
                curr_name = out.var(k).pretty_name;
                match = 1;
                break
            end
        end
        assert(match == 1)
    
        nonan = ~any([isnan(curr),isnan(curr_out)],2);
        curr_nonan = curr(nonan);
        curr_out_nonan = curr_out(nonan);
    
        if 0
            table(curr,curr_out)
        end
    
        
        tbl_2x2 = crosstab(curr_nonan,curr_out_nonan);
        [~,p,stats] = fishertest(tbl_2x2);
        
    
        if j == 1
            outcome_pred_text{i} = [outcome_pred_text{i},...
                {curr_name,sprintf('%1.2f (%1.2f-%1.2f)\nN=%d of %d with data',stats.OddsRatio,...
                stats.ConfidenceInterval(1),stats.ConfidenceInterval(2),sum(nonan),length(curr)),...
                sprintf('%1.2f',p)}];
        else
            outcome_pred_text{i} = [outcome_pred_text{i},...
                {sprintf('%1.2f (%1.2f-%1.2f)\nN=%d of %d with data',stats.OddsRatio,...
                stats.ConfidenceInterval(1),stats.ConfidenceInterval(2),sum(nonan),length(curr)),...
                sprintf('%1.2f',p)}];
        end
        
    end
end

combined_array = outcome_pred_text{1};
for i = 2:length(outcome_pred_text)
    combined_array = vertcat(combined_array, outcome_pred_text{i});
end

%{
table4T = cell2table(combined_array,'VariableNames',...
    {'Variable','OR (95% CI) for >50% event reduction (3 months)','p-value1',...
    'OR (95% CI) for event resolution (12-20 months)','p-value2',...
    'ED visits or hospitalizations','p-value3'});
%}
table4T = cell2table(combined_array,'VariableNames',...
    {'Variable','OR (95% CI) for event resolution (12-20 months)','p-value2',...
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
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',pre_median,pre_iqr(1),pre_iqr(2),sum(~isnan(pre_emu(group))),length(pre_emu(group))),...
            sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',post_median,post_iqr(1),post_iqr(2),sum(~isnan(post_emu(group))),length(post_emu(group))),...
            sprintf('T+: %1.1f',signedrank),...
            sprintf('%s',formatted_p_values(p))}];
        else
            table5_text{iphase} = [table5_text{iphase},{sprintf('%1.1f (%1.1f-%1.1f)\n(N = %d of %d with data)',post_median,post_iqr(1),post_iqr(2),...
                sum(~isnan(post_emu(group))),length(post_emu(group))),...
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
            all_supp_table_comp = [all_supp_table_comp;out.var(j).text_missing];
        end
    end
end
tableST1 = cell2table(all_supp_table_comp,'VariableNames',...
    {'Variable','All','Pre','Post','Statistic','p-value'});
writetable(tableST1,[results_folder,'TableS1.csv'],'QuoteStrings',true);

%% Prepare paired analyses
pre_freq_var = 'event_frequency_pre';
post_freq_vars = {'event_frequency_post','post_emu_freq_12mo_comb'};
npost = length(post_freq_vars);


%% Supp table 2 ED visits
tableST2_text = cell(2,1);
pre_ed_name = 'any_ed_pre';
post_ed_name = 'comb_ed_or_hospital_visits';

% Loop over two ways of defining post-EMU event frequency (visit vs survey)
for iphase = 1:2
    if iphase == 1
        group = pre;
        which_group = 'Care-as-usual';
    else
        group = post;
        which_group = 'Intervention';
    end

    % Add comparison for ED pre and post, using a Fisher's exact test to
    % compare the proportion of patients pre- and post-EMU with ED visits.
    pre_ed = T.(pre_ed_name)(group);
    post_ed = T.(post_ed_name)(group);

    % Remove any with nans
    anynan = any(isnan([pre_ed,post_ed]),2);
    pre_ed(anynan) = [];
    post_ed(anynan) = [];

    % do fisher exact test
    pre1 = sum(pre_ed==1);
    pre0 = sum(pre_ed==0);
    post1 = sum(post_ed == 1);
    post0 = sum(post_ed==0);
    tbl_2x2 = table([pre0;pre1],[post0;post1],...
        'VariableNames',{'Pre','Post'},'RowNames',{'0','1'});
    [~,p,stats] = fishertest(tbl_2x2);

    % prepare text
    
    tableST2_text{iphase} = [tableST2_text{iphase},{sprintf('%s',which_group),...
        sprintf('N = %d',sum(~isnan(pre_ed))),...
        sprintf('N = %d',sum(~isnan(post_ed))),...
        sprintf('OR (CI): %1.1f (%1.1f-%1.1f)',stats.OddsRatio,...
        stats.ConfidenceInterval(1),stats.ConfidenceInterval(2)),...
        sprintf('%s',formatted_p_values(p));...
        sprintf('No'),...
        sprintf('%d (%1.1f%%)',pre0,pre0/(pre0+pre1)*100),...
        sprintf('%d (%1.1f%%)',post0,post0/(post0+post1)*100),'','';...
        sprintf('Yes'),...
        sprintf('%d (%1.1f%%)',pre1,pre1/(pre0+pre1)*100),...
        sprintf('%d (%1.1f%%)',post1,post1/(post0+post1)*100),'',''}];
    
end

combined_array = tableST2_text{1};
for i = 2:length(tableST2_text)
    combined_array = vertcat(combined_array, tableST2_text{i});
end

table5T = cell2table(combined_array,'VariableNames',...
    {'Phase','Pre-EMU ED visits','Post-EMU ED visits',...
    'Odds-ratio','p-value'});
writetable(table5T,[results_folder,'tableS2.csv'],'QuoteStrings',true);

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

%% 3 point paired plot

figure
set(gcf,'position',[1 10 1440 500])
tiledlayout(1,2,'tilespacing','compact','padding','compact')
max_point = 40;

points = {'event_frequency_pre','event_frequency_post','post_emu_freq_12mo_comb'};

% Loop over arms
for ia = 1:2
    
    if ia == 1
        arm_text = 'Care-as-usual';
    else
        arm_text = 'Intervention';
    end
    nexttile
    
    % Get the data
    orig_pre = (T.(points{1})(T.phase==ia));
    orig_months3 = (T.(points{2})(T.phase==ia));
    orig_months12 = (T.(points{3})(T.phase==ia));
    npts = sum(T.phase==ia);

    % do the sign rank between pre and months3, and pre and months12
    p1 = signrank(orig_pre,orig_months3);
    p2 = signrank(orig_pre,orig_months12);

    %[pre,Ipre] = min([orig_pre,max_point]);
    %[months3,Imonths3] = min([orig_months3,max_point]);
    %[months12,Imonths12] = min([orig_months12,max_point]);

    pre = orig_pre;
    months3 = orig_months3;
    months12 = orig_months12;

    jitter = randn(npts,3)*0.03;
    
    fprintf('The following outlying points are omitted from the %s plot:\n',arm_text)

    % Loop over patients
    for ip = 1:npts
        all_points = [pre(ip) months3(ip) months12(ip)];
        if pre(ip) > max_point
            fprintf('Pre-EMU %d events/month\n',pre(ip));
        end
        if months3(ip) > max_point
            fprintf('3-months post %d events/month\n',months3(ip));
        end
        if months12(ip) > max_point
            fprintf('12-months post %d events/month\n',months12(ip));
        end

        % plot all individual points
        plot(1+jitter(ip,1),pre(ip),'ko','linewidth',2,'markersize',10)
        hold on
        plot(2+jitter(ip,2),months3(ip),'ko','linewidth',2,'markersize',10)
        plot(3+jitter(ip,3),months12(ip),'ko','linewidth',2,'markersize',10)

        % plot lines connecting the points
        % preferentially 1-2-3, but do 1-3 if there is no 2
        %{
        if ~isnan(months3(ip)) && ~isnan(months12(ip)) && ~isnan(pre(ip))
            plot([1+jitter(ip,1) 2+jitter(ip,2) 3+jitter(ip,3)],[pre(ip) months3(ip) months12(ip)],'k-','linewidth',2)
        elseif ~isnan(months3(ip)) && isnan(months12(ip)) && ~isnan(pre(ip))
            plot([1+jitter(ip,1) 2+jitter(ip,2)],[pre(ip) months3(ip)],'k-','linewidth',2)
        elseif isnan(months3(ip)) && ~isnan(months12(ip)) && ~isnan(pre(ip))
            plot([1+jitter(ip,1) 3+jitter(ip,3)],[pre(ip) months12(ip)],'k-','linewidth',2)
        elseif isnan(months3(ip)) && isnan(months12(ip)) && ~isnan(pre(ip))
        else
            error('why')
        end
        %}
        for k = 1:2
            if all_points(k) > max_point
                plot([k + jitter(ip,k) k+1+jitter(ip,k+1)],[max_point all_points(k+1)],'r--','linewidth',2)
            elseif all_points(k+1) > max_point
                plot([k + jitter(ip,k) k+1+jitter(ip,k+1)],[all_points(k) max_point],'r--','linewidth',2)
            else
                plot([k + jitter(ip,k) k+1+jitter(ip,k+1)],[all_points(k) all_points(k+1)],'k-','linewidth',2)
            end
        end

    end

    xticks([1 2 3])
    %{
    labels = cell(1,3);
    labels{1} = sprintf('Pre-EMU\n(N = %d)',sum(~isnan(pre)));
    labels{2} = sprintf('3 months post\n(N = %d)',sum(~isnan(months3)));
    labels{3} = sprintf('12-20 months post\n(N = %d)',sum(~isnan(months12)));
    %}
    %xticklabels({'Pre-EMU','3 months post','12-20 months post'})
    row1 = {'Pre-EMU','3 months post','12-20 months post'};
    row2 = {sprintf('(N = %d)',sum(~isnan(pre))),...
        sprintf('(N = %d)',sum(~isnan(months3))),...
        sprintf('(N = %d)',sum(~isnan(months12)))};
    labelArray = [row1; row2];
    %labelArray = strjust(pad(labelArray),'center');
    tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    xticklabels(tickLabels)
    ylabel('Events/month')
    title(arm_text)
    set(gca,'fontsize',20)
    
    ax.TickLabelInterpreter = 'tex';
    ylim([0 max_point])

    yl = ylim;
    ybar1 = yl(1) + (yl(2)-yl(1))* 1.05;
    ybar1_ticks = yl(1)+(yl(2)-yl(1))*1.03;
    ytext1 = yl(1) + (yl(2)-yl(1))* 1.09;
    ybar2 = yl(1) + (yl(2)-yl(1))* 1.14;
    ybar2_ticks = yl(1)+(yl(2)-yl(1))*1.12;
    ytext2 = yl(1) + (yl(2)-yl(1))* 1.18;
    yl_new = [yl(1) yl(1) + (yl(2)-yl(1))* 1.22];

    plot(xlim,[max_point, max_point],'k--')
    plot([1 2],[ybar1 ybar1],'k-','linewidth',2)
    plot([1 1],[ybar1_ticks ybar1],'k-','linewidth',2)
    plot([2 2],[ybar1_ticks ybar1],'k-','linewidth',2)
    ptext1 = def_ptext(p1);
    text(1.5,ytext1,ptext1,'fontsize',20,'HorizontalAlignment','center')

    plot([1 3],[ybar2 ybar2],'k-','linewidth',2)
    plot([1 1],[ybar2_ticks ybar2],'k-','linewidth',2)
    plot([3 3],[ybar2_ticks ybar2],'k-','linewidth',2)
    ptext2 = def_ptext(p2);
    text(2,ytext2,ptext2,'fontsize',20,'HorizontalAlignment','center')
    ylim(yl_new)
    yticks(0:5:40)


    
end
set(gcf,'renderer','painters')
print(gcf,[results_folder,'paired_plots3'],'-dpng')


