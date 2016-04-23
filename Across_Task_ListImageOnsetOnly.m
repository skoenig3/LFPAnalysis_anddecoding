%written 11/10 && 11/11 2015 by Seth Konig using/combining/modify code
%previously written to anlayze task speartely. Code focuses on eye movments
%or lack there off!!!!

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\'; %where to get data from
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Across Tasks Figures\'; %where to save figures

%sessions with recordings on the same day
listsq_files = {'PW140729_3-sorted.nex','PW140730_3-sorted.nex','PW140806_3-sorted.nex',...
    'PW140829_3-sorted.nex','PW141007_3-sorted.nex','PW141008_3-sorted.nex',...
    'PW141009_3-sorted.nex','PW141015_3-sorted.nex','PW141024_3-sorted.nex',...
    'PW141028_3-sorted.nex','PW150205_3-sorted.nex'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Important parametrs/values for ListSQ List
imgon_code = 23;
imgoff_code = 24;
fixspoton = 35;
imageX = 800; %horizontal size of the screen
imageY = 600; %horizontal size of the screen
minfixdur = 100;%how long of fixation is required after a saccade
maxoutside = 1000;%how long cumulatively have they looked outside the image before the data is ignored
%1000 gives them a reasonable buffer


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Important parameters/values across all tasks
trial_start_code = 15;
fixation_oncross = 8;
reward_code = 3;
Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity
Fs = 1000;
buffer1 = 2048;
buffer2 = 2048;


event_time_lock_data = cell(4,length(listsq_files));
low_freq_data = cell(4,length(listsq_files));
high_freq_data =  cell(4,length(listsq_files));
low_freq_coherence_data = cell(4,length(listsq_files));
high_freq_coherence_data = cell(4,length(listsq_files));

for file = 1:length(listsq_files)
    
    disp(['Running File #' num2str(file) ': ' listsq_files{file}(1:end-13) ])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---List portion of ListSQ---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import & Preprocess Data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([data_dir,listsq_files{file}(1:end-11)  '-preprocessed'],'cfg','fixationstats','item_set');
    
    %get important task specific information
    [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
    [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
    
    
    if str2double(item_set(7:8)) <= 6
        maximgdur = 7000*1.25;
        %images were up for 7 seconds cumulative
    else
        maximgdur = 5000*1.25;
        %images were up for 5 seconds cumulative
    end
    
    %preallocate space and parallel structure of cfg
    successful_trials = NaN(1,length(cfg.trl));
    for t = 1:length(cfg.trl);
        if any(cfg.trl(t).allval == imgon_code) %in which image was displayed or 1st item in sequence was displayed
            if itmlist(cfg.trl(t).cnd-1000) > sequence_items(end)% sequence trials and only want rewarded ones
                imgdur = cfg.trl(t).alltim(cfg.trl(t).allval == imgoff_code)-cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code);
                if imgdur < maximgdur %so didn't look away too too much, trying to only take the best data
                    successful_trials(t) = t;
                end
            end
        end
    end
    successful_trials = laundry(successful_trials);
    fixationstats = fixationstats(successful_trials);
    cfg.trl = cfg.trl(successful_trials);
    num_trials = length(successful_trials);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---process eye data and fixation durations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %only want trials in which we have a novel and repeat pair
    trials_with_img_pairs = NaN(1,num_trials);
    for t = 1:num_trials
        if sum(which_img == which_img(t)) == 2;
            trials_with_img_pairs(t) = t;
        end
    end
    trials_with_img_pairs = laundry(trials_with_img_pairs);
    
    l_trials(file) = length(trials_with_img_pairs)/2;
    if length(trials_with_img_pairs)/2 < 48
        disp(['Skipping File #' num2str(file) ': ' listsq_files{file}(1:end-11)])
        continue %need more data than this
    end
    
    fixationstats = fixationstats(trials_with_img_pairs);
    cfg.trl = cfg.trl(trials_with_img_pairs);
    which_img = which_img(trials_with_img_pairs);
    novel_vs_repeat = novel_vs_repeat(trials_with_img_pairs);
    num_trials = length(which_img);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---process eye data and fixation durations---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %only want trials in which we have a novel and repeat pair
    trials_with_img_pairs = NaN(1,num_trials);
    for t = 1:num_trials
        if sum(which_img == which_img(t)) == 2;
            trials_with_img_pairs(t) = t;
        end
    end
    trials_with_img_pairs = laundry(trials_with_img_pairs);
    
    if length(trials_with_img_pairs)/2 < 48
        continue %need more data than this
    end
    
    fixationstats = fixationstats(trials_with_img_pairs);
    cfg.trl = cfg.trl(trials_with_img_pairs);
    which_img = which_img(trials_with_img_pairs);
    novel_vs_repeat = novel_vs_repeat(trials_with_img_pairs);
    num_trials = length(which_img);
    
    all_event_times = NaN(num_trials,2);
    event_codes = [fixation_oncross imgon_code];
    for t = 1:num_trials
        trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
        for event = 1:length(event_codes);
            event_time = cfg.trl(t).alltim(cfg.trl(t).allval == event_codes(event));
            event_time = event_time(1);%fixation code can be found multiple times just want the first one
            all_event_times(t,event) = event_time-trial_start;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    
    nov_img_aligned = NaN(num_trials/2,3);
    rep_img_aligned = NaN(num_trials/2,3);
    nov_fix_aligned = NaN(num_trials/2,3);
    rep_fix_aligned = NaN(num_trials/2,3);
    
    nov_index = 1;
    rep_index = 1;
    for trl = 1:length(cfg.trl);
        startind = cfg.trl(t).begsmpind;%trial start index in all data
        if novel_vs_repeat(trl) == 1 %novel image
            
            nov_img_aligned(nov_index,:) = [all_event_times(t,1)-1-buffer1+startind ...
                all_event_times(t,2)-1+buffer2+startind -buffer1];
            nov_fix_aligned(nov_index,:) = [all_event_times(t,2)-1-buffer1+startind ...
                all_event_times(t,2)-1+buffer2+startind -buffer1];
            
            nov_index = nov_index+1;
        else %repeat image
            
            rep_img_aligned(rep_index,:) = [all_event_times(t,1)-1-buffer1+startind ...
                all_event_times(t,2)-1+buffer2+startind -buffer1];
            rep_fix_aligned(rep_index,:) = [all_event_times(t,2)-1-buffer1+startind ...
                all_event_times(t,2)-1+buffer2+startind -buffer1];
            
            rep_index = rep_index+1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Import LFP data---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %also filters out line noise
    
    %code to find spike/LFP channels
    LFPchannels = find_desired_channels(cfg,'LFP');
    
    % read the data from file and preprocess them
    cfg.channel       = cfg.channel(LFPchannels)';
    cfg.dftfilter     = 'yes';
    cfg.dftfreq       = Fline;
    cfg.padding       = 1;
    
    %LFP data aligned novel fixation on cross
    cfg.trl =  nov_img_aligned;
    novfix_aligneddata = ft_preprocessing(cfg);
    
    %LFP data aligned novel image onset
    cfg.trl = nov_fix_aligned;
    novimg_aligneddata = ft_preprocessing(cfg);
    
    %LFP data aligned repeat fixation on cross
    cfg.trl = rep_img_aligned;
    repfix_aligneddata = ft_preprocessing(cfg);
    
    %LFP data aligned repeat image onset
    cfg.trl = rep_fix_aligned;
    repimg_aligneddata = ft_preprocessing(cfg);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All Low Freq Power Spectrum---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfgfrq = [];
    cfgfrq.baseline = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(4,2);
    
    figure
    freq = lfp_powerspectrum(novfix_aligneddata,'low','all',[-0.75 -0.5]);
    low_freq_data{1,file} = freq;
    subplot(2,2,1)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Fixation On cross(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    freq = lfp_powerspectrum(repfix_aligneddata,'low','all',[-0.75 -0.5]);
    low_freq_data{2,file} = freq;
    subplot(2,2,2)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Fixation On cross(ms)')
    yl(1,:) = caxis;
    title('List Repeat Images')
    
    freq = lfp_powerspectrum(novimg_aligneddata,'low','all',[-1 -0.75]);
    low_freq_data{3,file} = freq;
    subplot(2,2,3)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Image Onset(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    freq = lfp_powerspectrum(repfix_aligneddata,'low','all',[-1 -0.75]);
    low_freq_data{4,file} = freq;
    subplot(2,2,4)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Image Onset(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:4
        subplot(2,2,sb)
        caxis([ymin ymax])
    end
    save_and_close_fig(figure_dir,['Across_Tasks_ListOnly_LowFreqPower' listsq_files{file}(1:end-13)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All high Freq Power Spectrum---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfgfrq = [];
    cfgfrq.baseline = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(4,2);
    
    figure
    freq = lfp_powerspectrum(novfix_aligneddata,'high','all',[-0.75 -0.5]);
    high_freq_data{1,file} = freq;
    subplot(2,2,1)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Fixation On cross(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    freq = lfp_powerspectrum(repfix_aligneddata,'high','all',[-0.75 -0.5]);
    high_freq_data{2,file} = freq;
    subplot(2,2,2)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Fixation On cross(ms)')
    yl(1,:) = caxis;
    title('List Repeat Images')
    
    freq = lfp_powerspectrum(novimg_aligneddata,'high','all',[-1 -0.75]);
    high_freq_data{3,file} = freq;
    subplot(2,2,3)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Image Onset(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    freq = lfp_powerspectrum(repfix_aligneddata,'high','all',[-1 -0.75]);
    high_freq_data{4,file} = freq;
    subplot(2,2,4)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Image Onset(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:4
        subplot(2,2,sb)
        caxis([ymin ymax])
    end
    save_and_close_fig(figure_dir,['Across_Tasks_ListOnly_highFreqPower' listsq_files{file}(1:end-13)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All high Freq Power Spectrum---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfgfrq = [];
    cfgfrq.baseline = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(4,2);
    
    figure
    freq = lfp_powerspectrum(novfix_aligneddata,'high','all',[-0.75 -0.5]);
    high_freq_data{1,file} = freq;
    subplot(2,2,1)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Fixation On cross(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    freq = lfp_powerspectrum(repfix_aligneddata,'high','all',[-0.75 -0.5]);
    high_freq_data{2,file} = freq;
    subplot(2,2,2)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Fixation On cross(ms)')
    yl(1,:) = caxis;
    title('List Repeat Images')
    
    freq = lfp_powerspectrum(novimg_aligneddata,'high','all',[-1 -0.75]);
    high_freq_data{3,file} = freq;
    subplot(2,2,3)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Image Onset(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    freq = lfp_powerspectrum(repfix_aligneddata,'high','all',[-1 -0.75]);
    high_freq_data{4,file} = freq;
    subplot(2,2,4)
    hold on
    ft_singleplotTFR(cfgfrq, freq);
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    hold off
    xlim([-0.75 0.75])
    xlabel('Time from Image Onset(ms)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for sb = 1:4
        subplot(2,2,sb)
        caxis([ymin ymax])
    end
    save_and_close_fig(figure_dir,['Across_Tasks_ListOnly_highFreqPower' listsq_files{file}(1:end-13)]);
    
end


