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
fixation_on_cross = 8;
reward_code = 3;
Fline = [59.8 59.9 60 60.1 60.2 119.8 119.9 120 120.1 120.2 179.8 179.9 180 180.1 180.2]; %line noise frequencies to remove
%60 Hz and it's harmonics as well as frequncies that are nearly equal to
%this as ther is some variabillity
Fs = 1000;
buffer1 = 1024;
buffer2 = 1024;


high_freq_coherence_data = cell(1,length(listsq_files));

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
                        if (fixationtimes(1,nextfix) - sacend) == 1 %should be == 1
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
    
    nov_saccade_aligned = NaN(4*num_trials,3);
    rep_saccade_aligned = NaN(4*num_trials,3);
    novfixation_duration = NaN(4*num_trials);
    repfixation_duration = NaN(4*num_trials);
    
    nov_index = 1;
    rep_index = 1;
    for trl = 1:length(cfg.trl);
        startind = cfg.trl(trl).begsmpind;%trial start index in all data
        
        sacs = saccade_times(trl,:);
        maxind = find(~isnan(sacs));
        for sac = 1:max(maxind)
            if ~isnan(saccade_times(trl,sac))
                if novel_vs_repeat(trl) == 1 %novel image
                    nov_saccade_aligned(nov_index,:) = [saccade_times(trl,sac)-1-buffer1+startind ...
                        saccade_times(trl,sac)-1+buffer2+startind -buffer1];
                    novfixation_duration(nov_index) = fixation_durations(trl,sac);
                    nov_index = nov_index+1;
                else %repeat image
                    rep_saccade_aligned(rep_index,:) = [saccade_times(trl,sac)-1-buffer1+startind ...
                        saccade_times(trl,sac)-1+buffer2+startind -buffer1];
                    repfixation_duration(nov_index) = fixation_durations(trl,sac);
                    rep_index = rep_index+1;
                end
            end
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
    
    %LFP data aligned saccades on novel images
    cfgnovsac = cfg;
    cfgnovsac.trl = nov_saccade_aligned;
    novsac_aligneddata = ft_preprocessing(cfgnovsac);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%---Plot All High Freq Coherence---%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    figure
    
    %Saccades during novel images during ListSQ List trials
    stat = lfp_phasecoherence2High(novsac_aligneddata,'all');
    high_freq_coherence_data{file} = stat;
    imagesc(stat.time,stat.freq,squeeze(abs(stat.cohspctrm(1,:,:))))
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    ylabel('Frequency (Hz)')
    xlim([-0.75 1.0])
    xlabel('Time from Saccade (sec)')
    yl(1,:) = caxis;
    title('List Novel Images')
    
    save_and_close_fig(figure_dir,['AcrossTask-highFreqHighTempresCoherence-' listsq_files{file}(1:end-13)]);
    
end
save('HighResTempCoherence_across_tasks')
emailme('Done with Across task data analysis')
%%
all_high_freq_coherence_data = [];
total_sessions = 0;
for file = 1:length(listsq_files)
    if ~isempty(high_freq_coherence_data{file})
        if file == 1
            all_high_freq_coherence_data = high_freq_coherence_data{file}.cohspctrm;
        else
            all_high_freq_coherence_data = all_high_freq_coherence_data+high_freq_coherence_data{file}.cohspctrm;
        end
        total_sessions = total_sessions + 1;
    end
end


figure

%Saccades during novel images during ListSQ List trials
stat = high_freq_coherence_data{1};
stat.cohspctrm = all_high_freq_coherence_data/total_sessions;
imagesc(stat.time,stat.freq,squeeze(abs(mean(stat.cohspctrm(:,:,:),1))))
line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
colorbar, axis xy
colormap jet
ylabel('Frequency (Hz)')
xlim([-0.75 1.0])
xlabel('Time from Saccade (sec)')

title('List Novel Images')