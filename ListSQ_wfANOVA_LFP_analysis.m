function ListSQ_wfANOVA_LFP_analysis(data_dir,preprocessed_data_file,figure_dir)
% written by Seth Konig August, 2014
% Function analyses LFPs from the ListSQ task using wavelet ANOVA analysis. 
% Function first removes like noise (60 Hz) and it's harmonics.
% Analysis is noparametric and makes no assumptions about the data!!!
%
% Inputs:
%   1) data_dir: direction where preprocessed data is located
%   2) preprocessed_data_file: preprocessed ListSQ data file with LFPs
%   divided by channel and trial. 
%   3) figure_dir: location of where to put save figures
%
% Outputs:
%   1) Saves figures to figure_dir 
%   2) Saves processed data to data_dir tagged with '-LFP_wfANOVA_results'

load([data_dir preprocessed_data_file],'data','hdr','cfg','item_set');
data_name = preprocessed_data_file(1:10);
[LFPchannels] = find_desired_channels(hdr,data,'LFP');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Remove line noise & it's harmonics---%%%
Fline       = 60; %removes harmonics so 120, 180, etc.
Fs = hdr.Fs; 

for trial = 1:length(data(LFPchannels(1)).values);
    LFPdata = [];
    for channel = 1:length(LFPchannels)
        LFPdata = [LFPdata; data(LFPchannels(channel)).values{trial}];
    end
    
    % taken directly from fieldtrip ft_preproc_dftfilter.m
    % determine the size of the data
    [Nchans, Nsamples] = size(LFPdata);
    
    % determine the largest integer number of line-noise cycles that fits in the data
    sel = 1:round(floor(Nsamples * Fline/Fs) * Fs/Fline);
    
    % temporarily remove mean to avoid leakage
    mdat = mean(LFPdata(:,sel),2);
    zeroedLFPdat  = LFPdata - mdat(:,ones(1,Nsamples));
    
    % fit a sin and cos to the signal and subtract them
    time  = (0:Nsamples-1)/Fs;
    tmp  = exp(j*2*pi*Fline*time);                    % complex sin and cos
    % ampl = 2*dat*tmp'/Nsamples;                  % estimated amplitude of complex sin and cos
    ampl = 2*zeroedLFPdat(:,sel)*tmp(sel)'/length(sel);     % estimated amplitude of complex sin and cos on integer number of cycles
    est  = ampl*tmp;                               % estimated signal at this frequency
    %filt = dat - est;                              % subtract estimated signal
    LFPdata = LFPdata - est; %do want to remove the mean
    LFPdata = real(LFPdata);
    LFPdata = LFPdata./(std(LFPdata')'*ones(1,length(LFPdata))); %normalize since impedance changes absolute power
    
    for channel = 1:4;
        data(LFPchannels(channel)).values{trial} = LFPdata(channel,:);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Parse LFPs by sequence and by novel vs repeat images---%%%
 [itmlist,sequence_items] = read_ListSQ_itm_and_cnd_files(item_set);

crosses_on_off = [23 25 27 29;
    24 26 28 30];
reward_code = 3;
img_on_code =23;

parsed_SQLFP_data = NaN(424*4*4,512);%trial * item in sequence*LFP channel
parsed_imgLFP_data = NaN(96*2*4,2048);%trial * item in sequence*LFP channel
seq = NaN(1,424*4*4);
trial_num = NaN(1,424*4*4);
which_channel = NaN(1,424*4*4);
item_num = NaN(1,424*4*4); %item in sequence
which_img = NaN(1,96*2*4);
which_img_channel = NaN(1,96*2*4);
trial_count = ones(1,2);
img_trial_num = 1;
img_trial_count = 1;
row = 1;
img_row = 1;
for t = 1:length(data(LFPchannels(1)).values);
    if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(2) %sequence trials
        if  any(cfg.trl(t).allval == reward_code);
            eyedatstart = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
            tests_on_off = NaN(2,4);
            for c = 1:size(crosses_on_off,2);
                tests_on_off(1,c) = cfg.trl(t).alltim(cfg.trl(t).allval == crosses_on_off(1,c))-eyedatstart;
                tests_on_off(2,c) = cfg.trl(t).alltim(cfg.trl(t).allval == crosses_on_off(2,c))-eyedatstart;
            end
            for sequence_item = 1:4
                for channel = 1:4;
                    parsed_SQLFP_data(row,:) = data(LFPchannels(channel)).values{t}...
                        (tests_on_off(2,sequence_item)-511:tests_on_off(2,sequence_item));
                    which_channel(row) = channel;
                    item_num(row) = sequence_item;
                    if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(1) %image trials ignore crosshair
                        trial_num(row)= trial_count(1);
                        seq(row) = 1;
                    else
                        trial_num(row)= trial_count(2);
                        seq(row) = 2;
                    end
                    row = row+1;
                end
            end
            if itmlist(cfg.trl(t).cnd-1000) <= sequence_items(1) %image trials ignore crosshair
                trial_count(1) = trial_count(1) + 1;
            else
                trial_count(2) = trial_count(2) + 1;
            end
            
        end
    else
        if  any(cfg.trl(t).allval == img_on_code);
            for channel = 1:4;
                img_ind_on = cfg.trl(t).alltim(cfg.trl(t).allval == img_on_code)-cfg.trl(t).alltim(cfg.trl(t).allval == 15);
                parsed_imgLFP_data(img_row,:) = data(LFPchannels(channel)).values{t}(img_ind_on:img_ind_on+2047);
                which_img(img_row) = itmlist(cfg.trl(t).cnd-1000);
                which_img_channel(img_row) = channel;
                img_row = img_row+1;
            end
            img_trial_count = img_trial_count+1;
        end
    end
end

% remove extra NaNs for Sequence portion
nanind = isnan(seq);
seq(nanind) = [];
trial_num(nanind) = [];
which_channel(nanind) = [];
item_num(nanind) = [];
parsed_SQLFP_data(nanind,:) = [];

novel_vs_repeat = NaN(1,96*2*4);
for img = 1:max(which_img)
    imgind = find(which_img == img);
    dimgind = find(diff(imgind) > 1);
    if ~isempty(imgind);
        novel_vs_repeat(imgind(1:dimgind)) = 1;
        novel_vs_repeat(imgind(dimgind+1:end)) = 2;
    end
end

% remove extra NaNs for Image portion
nanind = isnan(which_img);
which_img(nanind) = [];
which_img_channel(nanind) = [];
novel_vs_repeat(nanind) = [];
parsed_imgLFP_data(nanind,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Run various wfANOVA analysis on LFPs---%%%

% Determine if there is a difference between sequences
sequence_contrasts = LFP_wavelet_ANOVA(parsed_SQLFP_data,seq);
save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference between Sequences'])

% Determine if there is a difference between sequences ignoring 1st item
tempdata = parsed_SQLFP_data;
tempdata(item_num == 1,:) = [];
tempseq = seq;
tempseq(item_num == 1) = [];
sequence_contrasts_minus_1st_item = LFP_wavelet_ANOVA(tempdata,tempseq);
save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference between Sequences ignoring item 1'])

% Determine if the there is a difference between items irrespective of sequence
item_contrast = LFP_wavelet_ANOVA(parsed_SQLFP_data,item_num,'Item');
save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference between Items irrespective of sequence'])

% Determine if there is a difference between sequences by item
item_by_sequence_contrasts = cell(1,4);
for item = 1:4;
    tempdata = parsed_SQLFP_data;
    tempdata = tempdata(item_num == item,:);
    tempseq = seq(item_num == item);
    item_by_sequence_contrasts{item}=LFP_wavelet_ANOVA(tempdata,tempseq);
    save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference in LFP in Sequence 1 vs Sequence 2 for Item ' num2str(item)])
end

% Determine if the there is a difference between items in sequence 1 and sequence 2
items_in_sequence_1_contrasts = LFP_wavelet_ANOVA(parsed_SQLFP_data(seq==1,:),item_num(seq==1),'Item');
save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference between Items in sequence 1'])

items_in_sequence_2_contrasts = LFP_wavelet_ANOVA(parsed_SQLFP_data(seq==2,:),item_num(seq==2),'Item');
save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference between Items in sequence 2'])

% Determine if there is a difference between recording electrodes/channels
channel_contrast = LFP_wavelet_ANOVA(parsed_SQLFP_data,which_channel);
save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference between Channels'])

% Determine if LFPs reflect novel vs repeat images
novel_vs_repeat_contrasts = LFP_wavelet_ANOVA(parsed_imgLFP_data,novel_vs_repeat);
save_and_close_fig(figure_dir,[data_name '_LFPwfANOVA_Difference between Novel and Repeat Images'])

save([data_dir data_name(1:10)  '-LFP_wfANOVA_results.m'],'sequence_contrasts','sequence_contrasts_minus_1st_item',...
    'item_contrast','item_by_sequence_contrasts','items_in_sequence_1_contrasts',....
    'items_in_sequence_2_contrasts','channel_contrast','novel_vs_repeat_contrasts');

end