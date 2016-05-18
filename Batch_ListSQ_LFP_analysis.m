% written by Seth Konig August 2014. Updated to V2 by SDK 1/7/16. Updated
% to import files, mutliunit stats, etc. from excel file. V2 also analyzes
% firing rates over time and uses only stable portion of the task. 
% code runs all the other code preprocess then process recording data
clar 

%monkey = 'Vivian'; 
monkey = 'Tobii'; 

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Figures\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Read in Excel Sheet for Session data---%%%
%only need to run when somethings changed or sessions have been added
if strcmpi(monkey,'Vivian')
    excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
    excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
    eye_data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
  
    load([eye_data_dir 'Across_Session_Unit_Data_Vivian.mat'])
    
    predict_rt = 156;%156 ms prediction 5-percentile
    chamber_zero = [13.5 -11]; %AP ML

elseif strcmpi(monkey,'Tobii')
    excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
    excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
    eye_data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';

    predict_rt = 138;%ms prediction 5-percentile
    chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
   
    load([eye_data_dir 'Across_Session_Unit_Data_Tobii.mat'])
    session_data(end-1) = [];%last file doesn't have strobe signal working on importing the data
end

% for sess = length(session_data)
%     get_saccade_aligned_LFPs_across_tasks(session_data{sess},data_dir,eye_data_dir)
% end

for sess = length(session_data)
    saccade_aligned_LFP_analysis(session_data{sess},data_dir,figure_dir)
end

emailme('Done Importing and processsing Tobiis LFP data')
%%
clar 

monkey = 'Vivian'; 
%monkey = 'Tobii'; 

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Figures\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---Read in Excel Sheet for Session data---%%%
%only need to run when somethings changed or sessions have been added
if strcmpi(monkey,'Vivian')
    excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Vivian\';
    excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
    eye_data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
  
    load([eye_data_dir 'Across_Session_Unit_Data_Vivian.mat'])
    
    predict_rt = 156;%156 ms prediction 5-percentile
    chamber_zero = [13.5 -11]; %AP ML

elseif strcmpi(monkey,'Tobii')
    excel_dir = '\\research.wanprc.org\Research\Buffalo Lab\eblab\PLX files\Tobii\';
    excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
    eye_data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';

    predict_rt = 138;%ms prediction 5-percentile
    chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
   
    load([eye_data_dir 'Across_Session_Unit_Data_Tobii.mat'])
    session_data(end) = [];%last file doesn't have strobe signal working on importing the data
end

for sess = 1:length(session_data)
    get_saccade_aligned_LFPs_across_tasks(session_data{sess},data_dir,eye_data_dir)
end

for sess = 1:length(session_data)
    saccade_aligned_LFP_analysis(session_data{sess},data_dir,figure_dir)
end
emailme('Done Importing and processsing Vivians LFP data')

