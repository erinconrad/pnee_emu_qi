%% ReadMe
%{
To run this code, do the following:

1) Download or clone the codebase: https://github.com/erinconrad/pnee_emu_qi
2) Download myFisher23 https://www.mathworks.com/matlabcentral/fileexchange/15399-myfisher23
and add it to your path (this is needed to do a Fisher exact test for
gender, which is a 2x3 table)
3) Create a Matlab script file called pnee_locations.m and put it in your
path. This should point to the data folder and the results folder. Here is
an example of how this should look:

    function locations = pnee_locations

    locations.main_folder = '/***path containing the codebase***/';
    locations.data = [locations.main_folder,'data/'];
    locations.scripts = [locations.main_folder,'scripts/'];
    locations.results = [locations.main_folder,'results/'];

    end
4) Download the erin_analysis_full csv file from redcap and put it in your
data folder (***eventually put this on PennBox so that anyone can run this)
5) Navigate to the script folder and type
    >> main

This will run the analyses and output to the results folder:
- pre_vs_post_table.csv, a table with the pre- vs post-intervention arm
comparisons
- paired_table.csv, a table with the paired pre-EMU vs post-EMU sz
frequency comparisons
- paired_plots.png, a figure with the paired pre-EMU vs post-EMU sz
frequency comparisons


%}