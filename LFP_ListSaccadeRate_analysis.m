%written 11/10 && 11/11 2015 by Seth Konig using/combining/modify code
%previously written to anlayze task speartely. Code focuses on eye movments
%or lack there off!!!!

clar
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\'; %where to get data from
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\SaccadeRate\'; %where to save figures

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
fixation_on_cross = 8;
reward_code = 3;
Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity
Fs = 1000;
buffer1 = 2048;
buffer2 = 2048;


fixdurthresh = [100 200;
    200 300;
    300 400;
    400 500
    500  750];

event_time_lock_data = cell(size(fixdurthresh,1),length(listsq_files));
low_freq_data = cell(size(fixdurthresh,1),length(listsq_files));
high_freq_data = cell(size(fixdurthresh,1),length(listsq_files));
low_freq_coherence_data = cell(size(fixdurthresh,1),length(listsq_files));
high_freq_coherenc_data = cell(size(fixdurthresh,1),length(listsq_files));
total_trials_per_session = cell(1,size(fixdurthresh,1));



for file = 1:length(listsq_files)
    
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
    
    saccade_times = NaN(num_trials,75);
    fixation_durations = NaN(num_trials,75); %keep track of fixation duration following saccade because we may want to sort things out
    for t = 1:length(fixationstats);
        if any(cfg.trl(t).allval == imgon_code) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
            trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code);
            imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == imgon_code);%to ignore 1st fixation
            imgoff =  cfg.trl(t).alltim(cfg.trl(t).allval == imgoff_code);%to ignore last fixation if cutoff
            saccadetimes = fixationstats{t}.saccadetimes;
            fixationtimes = fixationstats{t}.fixationtimes;
            if ~isempty(fixationtimes)
                %ignore saccades before the first 250 ms of image on to
                %remove image onset confounds
                
                if saccadetimes(1,1) < fixationtimes(1,1) %saccade is the first index
                    sacindex = 0;
                else
                    sacindex = +1;
                end
                timeoutside = cumsum(isnan(fixationstats{t}.XY(1,:)));
                timeoutside = find(timeoutside >= maxoutside); %looked away a lot, again trying to cut out more data but with reason
                if ~isempty(timeoutside)
                    timeoutside = timeoutside(1);
                    toolate = find(saccadetimes(1,:) > timeoutside);
                    saccadetimes(:,toolate) = NaN;
                    toolate = find(fixationtimes(1,:) > timeoutside);
                    fixationtimes(:,toolate) = NaN;
                end
                toolate = find(saccadetimes(2,:)) > imgoff-trial_start;
                saccadetimes(:,toolate) = NaN;
                toolate = find(fixationtimes(2,:)) > imgoff-trial_start;
                fiationtimes(:,toolate) = NaN;
                
                tooearly = find(saccadetimes(1,:) <= imgon-trial_start+250);
                saccadetimes(:,tooearly) = NaN;
                tooearly = find(fixationtimes(1,:) <= imgon-trial_start+250);
                fixationtimes(:,tooearly) = NaN;
                
                %time outside image is ignored so not 1:1 saccaddes to
                %fixations
                fixdur =  NaN(1,size(saccadetimes,2));
                for sac = 1:size(saccadetimes,2);
                    sacend = saccadetimes(2,sac);
                    nextfix = find(fixationtimes(1,:) > sacend);
                    if ~isempty(nextfix)
                        nextfix = nextfix(1);
                        if (fixationtimes(1,nextfix) - sacend) == 1 %should be == 1 unless saccade is outside of image
                            fixdur(sac) = fixationtimes(2,nextfix)-fixationtimes(1,nextfix)+1;
                        end
                    end
                end
                
                goodsacs = find(~isnan(saccadetimes(1,:)) & fixdur > minfixdur);
                saccade_times(t,1:length(goodsacs)) = saccadetimes(1,goodsacs);
                fixation_durations(t,1:length(goodsacs)) = fixdur(goodsacs);
            end
        end
    end
    sumnan = nansum(saccade_times);
    badcols = find(sumnan == 0);
    saccade_times(:,badcols)  = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---convert Data into Field Trip Friendly Format---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % want following organization [event start index-buffer1, event end index+buffer2, -buffer1]
    
    novtrials = find(novel_vs_repeat == 1);
    nov_saccade_aligned = cell(1,size(fixdurthresh,1));
    novfixation_duration = cell(1,size(fixdurthresh,1));
    for fdt = 1:size(fixdurthresh,1)
        total_fix = sum(sum(fixation_durations(novtrials,:) >= fixdurthresh(fdt,1) & fixation_durations(novtrials,:) < fixdurthresh(fdt,2)));
        nov_saccade_aligned{fdt} = NaN(total_fix,3);
        novfixation_duration{fdt} = NaN(1,total_fix);
    end
    
    
    nov_index = ones(1,size(fixdurthresh,1));
    for trl = 1:length(cfg.trl);
        startind = cfg.trl(trl).begsmpind;%trial start index in all data
        
        sacs = saccade_times(trl,:);
        maxind = find(~isnan(sacs));
        for sac = 1:max(maxind)
            if ~isnan(saccade_times(trl,sac))
                if novel_vs_repeat(trl) == 1 %novel image
                    which_index = find(fixation_durations(trl,sac) >= fixdurthresh(:,1) & fixation_durations(trl,sac) < fixdurthresh(:,2));
                    if ~isempty(which_index) %should only happen if fixation duration is greater than max
                        nov_saccade_aligned{which_index}(nov_index(which_index),:) = [saccade_times(trl,sac)-1-buffer1+startind ...
                            saccade_times(trl,sac)-1+buffer2+startind -buffer1];
                        novfixation_duration{which_index}(nov_index(which_index)) = fixation_durations(trl,sac);
                        nov_index(which_index) = nov_index(which_index)+1;
                    end
                end
            end
        end
    end
    total_trials_per_session{file} =  nov_index-1;
    
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
    
    %LFP data aligned to saccades during novel images
    
    novsac_aligneddata = cell(1,size(fixdurthresh,1));
    cfgnovsac = cell(1,size(fixdurthresh,1));
    for ftd = 1:size(fixdurthresh,1);
        cfgnovsac{ftd} = cfg;
        cfgnovsac{ftd}.trl = nov_saccade_aligned{ftd};
        novsac_aligneddata{ftd} = ft_preprocessing(cfgnovsac{ftd});
    end
    
    %---Plot Average Saccade Locked LFP---%
    yl = NaN(size(fixdurthresh,1),2);
    figure
    for ftd = 1:size(fixdurthresh,1);
        timelock = ft_timelockanalysis(cfgnovsac{ftd},novsac_aligneddata{ftd});
        event_time_lock_data{ftd,file} = timelock;
        subplot(2,3,ftd)
        plot(timelock.time,timelock.avg')
        xlim([-0.5 1.5])
        xlabel('Time from Saccade (sec)')
        ylabel('Voltage (uV)')
        yl(ftd,:) = ylim;
        title(['Fixation durations: ' num2str(fixdurthresh(ftd,1)) '-' num2str(fixdurthresh(ftd,2))...
            ' ms (n =' num2str(nov_index(ftd)-1) ')'])
    end
    
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for ftd = 1:size(fixdurthresh,1);
        subplot(2,3,ftd)
        hold on
        plot([0 0],[ymin ymax],'k--')
        hold off
        set(gca,'XMinorTick','on','YMinorTick','on')
        grid on
        grid(gca,'minor')
        ylim([ymin ymax])
    end
    save_and_close_fig(figure_dir,['SaccadeRate_LFP_' listsq_files{file}(1:end-11)])
    
    %---Plot Low Frequency Power Spectrum---%
    cfgfrq = [];
    cfgfrq.baseline = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(6,2);
    figure
    for ftd = 1:size(fixdurthresh,1);
        %Saccades during novel images during ListSQ List trials
        freq = lfp_powerspectrum(novsac_aligneddata{ftd},'low','all',[-0.75 -0.5]);
        low_freq_data{ftd,file} = freq;
        subplot(2,3,ftd)
        ft_singleplotTFR(cfgfrq,freq);
        xlim([-0.5 1.5])
        xlabel('Time from Saccade (sec)')
        yl(ftd,:) = caxis;
        title(['Fixation durations: ' num2str(fixdurthresh(ftd,1)) '-' num2str(fixdurthresh(ftd,2))...
            ' ms (n =' num2str(nov_index(ftd)-1) ')'])
    end
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for ftd = 1:size(fixdurthresh,1);
        subplot(2,3,ftd)
        hold on
        plot([0 0],[3.5 30.5],'w-')
        hold off
        caxis([ymin ymax])
    end
    save_and_close_fig(figure_dir,['SaccadeRate_LowFreqPowerSpectrum_' listsq_files{file}(1:end-11)])
    
    %---Plot High Frequency Power Spectrum---%
    cfgfrq = [];
    cfgfrq.baseline = [-0.75 -0.5];
    cfgfrq.baselinetype = 'absolute';
    cfgfrq.maskstyle    = 'saturation';
    cfgfrq.zparam       = 'powspctrm';
    
    yl = NaN(6,2);
    figure
    for ftd = 1:size(fixdurthresh,1);
        %Saccades during novel images during ListSQ List trials
        freq = lfp_powerspectrum(novsac_aligneddata{ftd},'high','all',[-0.75 -0.5]);
        high_freq_data{ftd,file} = freq;
        subplot(2,3,ftd)
        ft_singleplotTFR(cfgfrq,freq);
        xlim([-0.5 1.5])
        xlabel('Time from Saccade (sec)')
        yl(ftd,:) = caxis;
        title(['Fixation durations: ' num2str(fixdurthresh(ftd,1)) '-' num2str(fixdurthresh(ftd,2))...
            ' ms (n =' num2str(nov_index(ftd)-1) ')'])
    end
    
    %rescale and add line plot
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for ftd = 1:size(fixdurthresh,1);
        subplot(2,3,ftd)
        hold on
        plot([0 0],[29 181],'w-')
        hold off
        caxis([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['SaccadeRate_HighFreqPowerSpectrum_' listsq_files{file}(1:end-11)])
    
   
    %---Plot Low Frequency Coherence---%
    yl = NaN(6,2);
    figure
    for ftd = 1:size(fixdurthresh,1);
        %Saccades during novel images during ListSQ List trials
        stat = lfp_phasecoherence(novsac_aligneddata{ftd},'all');
        low_freq_coherence_data{ftd,file} = stat;
        subplot(2,3,ftd)
        imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
        line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
        colorbar, axis xy
        colormap jet
        ylabel('Frequency (Hz)')
        xlim([-0.5 1.5])
        xlabel('Time from Saccade (sec)')
        yl(ftd,:) = caxis;
        title(['Fixation durations: ' num2str(fixdurthresh(ftd,1)) '-' num2str(fixdurthresh(ftd,2))...
            ' ms (n =' num2str(nov_index(ftd)-1) ')'])
    end
    
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for ftd = 1:size(fixdurthresh,1);
        subplot(2,3,ftd)
        caxis([ymin ymax])
    end
    
    save_and_close_fig(figure_dir,['SaccadeRate_LowFreqCoherence_' listsq_files{file}(1:end-11)])
    
    
    %---Plot High Frequency Coherence---%
    yl = NaN(6,2);
    figure
    for ftd = 1:size(fixdurthresh,1);
        %Saccades during novel images during ListSQ List trials
        stat = lfp_phasecoherence2High(novsac_aligneddata{ftd},'all');
        high_freq_coherence_data{ftd,file} = stat;
        subplot(2,3,ftd)
        imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
        line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
        colorbar, axis xy
        colormap jet
        ylabel('Frequency (Hz)')
        xlim([-0.5 1.5])
        xlabel('Time from Saccade (sec)')
        yl(ftd,:) = caxis;
        title(['Fixation durations: ' num2str(fixdurthresh(ftd,1)) '-' num2str(fixdurthresh(ftd,2))...
            ' ms (n =' num2str(nov_index(ftd)-1) ')'])
    end
    
    ymin = min(yl(:,1));
    ymax = max(yl(:,2));
    for ftd = 1:size(fixdurthresh,1);
        subplot(2,3,ftd)
        caxis([ymin ymax])
    end
    save_and_close_fig(figure_dir,['SaccadeRate_HighFreqCoherence_' listsq_files{file}(1:end-11)])

end
save('Saccade_Rate_Data')
emailme('Saccade Rate code is done')
%%
total_trials_per_session{11} = []; %just for filler

session_avg_event_time_lock_data = cell(1,size(fixdurthresh,1));
session_avg_low_freq_data = cell(1,size(fixdurthresh,1));
session_avg_high_freq_data = cell(1,size(fixdurthresh,1));
session_avg_low_freq_coherence_data = cell(1,size(fixdurthresh,1));
session_avg_high_freq_coherenc_data = cell(1,size(fixdurthresh,1));

session_total_samples = zeros(1,size(fixdurthresh,1));
total_sessions = 0;
for file = 1:length(listsq_files)
    for fdt = 1:size(fixdurthresh,1)
        if ~isempty(event_time_lock_data{fdt,file})
            if file == 1
                session_avg_event_time_lock_data{fdt} = event_time_lock_data{fdt,file}.avg;
                session_avg_low_freq_data{fdt} = low_freq_data{fdt,file}.powspctrm;
                session_avg_high_freq_data{fdt} = high_freq_data{fdt,file}.powspctrm;
                session_avg_low_freq_coherence_data{fdt} = low_freq_coherence_data{fdt,file}.cohspctrm;
                session_avg_high_freq_coherence_data{fdt} = high_freq_coherence_data{fdt,file}.cohspctrm;
            else
                session_avg_event_time_lock_data{fdt} = session_avg_event_time_lock_data{fdt}+event_time_lock_data{fdt,file}.avg;
                session_avg_low_freq_data{fdt} = session_avg_low_freq_data{fdt}+low_freq_data{fdt,file}.powspctrm;
                session_avg_high_freq_data{fdt} = session_avg_high_freq_data{fdt}+high_freq_data{fdt,file}.powspctrm;
                session_avg_low_freq_coherence_data{fdt} = session_avg_low_freq_coherence_data{fdt}+low_freq_coherence_data{fdt,file}.cohspctrm;
                session_avg_high_freq_coherence_data{fdt} = session_avg_high_freq_coherence_data{fdt}+high_freq_coherence_data{fdt,file}.cohspctrm;
            end
        end
    end
    if ~isempty(total_trials_per_session{file})
        session_total_samples = session_total_samples+total_trials_per_session{file};
        total_sessions = total_sessions+1;
    end
end


%---Plot time locked averge LFP across sessions---%

yl = NaN(size(fixdurthresh,1),2);
figure
for fdt = 1:size(fixdurthresh,1);
    subplot(2,3,fdt)
    plot(event_time_lock_data{1}.time,session_avg_event_time_lock_data{fdt}'/total_sessions)
    xlim([-0.5 1.5])
    xlabel('Time from Saccade (sec)')
    ylabel('Voltage (uV)')
    yl(fdt,:) = ylim;
    title(['Fixation durations: ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2))...
        ' ms (n =' num2str(session_total_samples(fdt)-1) ')'])
end

ymin = 0.75*min(yl(:,1));
ymax = 0.75*max(yl(:,2));
for fdt = 1:size(fixdurthresh,1);
    subplot(2,3,fdt)
    hold on
    plot([0 0],[ymin ymax],'k--')
    hold off
    set(gca,'XMinorTick','on','YMinorTick','on')
    grid on
    grid(gca,'minor')
    ylim([ymin ymax])
end
save_and_close_fig(figure_dir,['Session_SaccadeRate_LFP_' listsq_files{file}(1:end-11)])


%---Plot Low Frequency Power Spectrum---%
cfgfrq = [];
cfgfrq.baseline = [-0.75 -0.5];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

yl = NaN(6,2);
figure
for fdt = 1:size(fixdurthresh,1);
    subplot(2,3,fdt)
    freq = low_freq_data{1};
    freq.powspctrm = session_avg_low_freq_data{fdt}/total_sessions;
    ft_singleplotTFR(cfgfrq,freq);
    xlim([-0.5 1.5])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(['Fixation durations: ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2))...
        ' ms (n =' num2str(session_total_samples(fdt)-1) ')'])
end

%rescale and add line plot
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for fdt = 1:size(fixdurthresh,1);
    subplot(2,3,fdt)
    hold on
    plot([0 0],[3.5 30.5],'w-')
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_SaccadeRate_LowFreqPowerSpectrum_' listsq_files{file}(1:end-11)])

%---Plot High Frequency Power Spectrum---%
cfgfrq = [];
cfgfrq.baseline = [-0.75 -0.5];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';

yl = NaN(6,2);
figure
for fdt = 1:size(fixdurthresh,1);
    %Saccades during novel images during ListSQ List trials
    freq = high_freq_data{1};
    freq.powspctrm = session_avg_high_freq_data{fdt}/total_sessions;
    subplot(2,3,fdt)
    ft_singleplotTFR(cfgfrq,freq);
    xlim([-0.5 1.5])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(['Fixation durations: ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2))...
        ' ms (n =' num2str(session_total_samples(fdt)-1) ')'])
end

%rescale and add line plot
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for fdt = 1:size(fixdurthresh,1);
    subplot(2,3,fdt)
    hold on
    plot([0 0],[29 181],'w-')
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_SaccadeRate_HighFreqPowerSpectrum_' listsq_files{file}(1:end-11)])


%---Plot Low Frequency Coherence---%
yl = NaN(6,2);
figure
for fdt = 1:size(fixdurthresh,1);
    %Saccades during novel images during ListSQ List trials
    stat = low_freq_coherence_data{1};
    stat.cohspctrm = session_avg_low_freq_coherence_data{fdt}/total_sessions;
    subplot(2,3,fdt)
    imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.5 1.5])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(['Fixation durations: ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2))...
        ' ms (n =' num2str(session_total_samples(fdt)-1) ')'])
end

ymin = min(yl(:,1));
ymax = max(yl(:,2));
for fdt = 1:size(fixdurthresh,1);
    subplot(2,3,fdt)
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_SaccadeRate_LowFreqCoherence_' listsq_files{file}(1:end-11)])

%---Plot High Frequency Coherence---%
yl = NaN(6,2);
figure
for fdt = 1:size(fixdurthresh,1);
    %Saccades during novel images during ListSQ List trials
    stat = high_freq_coherence_data{1};
    stat.cohspctrm = session_avg_high_freq_coherence_data{fdt}/total_sessions;
    subplot(2,3,fdt)
    imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.5 1.5])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(['Fixation durations: ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2))...
        ' ms (n =' num2str(session_total_samples(fdt)-1) ')'])
end

ymin = min(yl(:,1));
ymax = max(yl(:,2));
for fdt = 1:size(fixdurthresh,1);
    subplot(2,3,fdt)
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,['Session_SaccadeRate_HighFreqCoherence_' listsq_files{file}(1:end-11)])
