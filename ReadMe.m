%% ReadMe
%{
To run this code, do the following:

1) Download or clone the codebase: https://github.com/erinconrad/pnee_emu_qi
2) Download myFisher23 https://www.mathworks.com/matlabcentral/fileexchange/15399-myfisher23
and add it to your path (this is needed to do a Fisher exact test for
gender, which is non-binary and so requires a 2x3 table)
3) Create a Matlab script file called pnee_locations.m and put it in your
path. This should point to the data folder and the results folder. Here is
an example of how this should look:

    function locations = pnee_locations

    locations.main_folder = '/***path containing the codebase***/';
    locations.data = [locations.main_folder,'data/'];
    locations.scripts = [locations.main_folder,'scripts/'];
    locations.results = [locations.main_folder,'results/'];

    end
4) Download the pnee_data.csv and pnee_labels.csv files from the following
link and add them to your data folder: 

https://upenn.box.com/s/h28wvjlj4381z3ddgwhkp0avk8iepcsn

Note for Erin: this data is from the redcap report called
erin_this_is_the_one_for_the_paper in PNEE EMU Pathway redcap

5) Navigate to the script folder and type
    >> main

This will run the analyses and output tables and the figure to the results folder.


Kelly Boylan and Erin Conrad
July 2024

%}