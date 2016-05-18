% % written Seth Konig 4/29/16
% % grab analyzed data across all sessions and monkeys and get the averages
% % across all sessions and do some statistical analysis
%
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\';

%fixation durations by category minimum already set at 100 ms
fixdurthresh = [100 200;
    200 300;
    300 400;
    400 500
    500  750];
Fs = 1000; %sanmpling frequence
smval = []; %no smoothing but want dofill plots
time_window = -0.75:.01:1.5; %what window to consider eye data over
baseline_time_interval = [0.25 0.75]; %in ITI

a = what(data_dir);
a = a.mat;

files = [];
for f = 1:size(a,1)
    if ~isempty(strfind(a{f},'preproccessed_saccade_aligned_LFPs'))
        files = [f files];
    end
end
%%
% %---Preallocate structure for variables---%
% list_average_LFP = cell(size(fixdurthresh,1),length(files));
% list_low_freq = cell(size(fixdurthresh,1),length(files));
% list_high_freq = cell(size(fixdurthresh,1),length(files));
% list_low_freq_coh = cell(size(fixdurthresh,1),length(files));
% list_high_freq_coh = cell(size(fixdurthresh,1),length(files));
%
%
% sequence_average_LFP = cell(1,length(files));
% sequence_low_freq = cell(1,length(files));
% sequence_high_freq = cell(1,length(files));
% sequence_low_freq_coh = cell(1,length(files));
% sequence_high_freq_coh = cell(1,length(files));
%
% sequence_average_LFP = cell(1,length(files));
% for f = 1:length(files)
%     disp([ 'Processing file ' a{files(f)}(1:8)])
%     load([data_dir a{files(f)}],'list_saccade_aligneddata','stretched_fixation_durations',...
%         'sequence_saccade_aligneddata','stretched_time_to_fixation','listsq_ITI_aligneddata');
%
%     %in case too_late was implemented too late
%     stretched_fixation_durations = stretched_fixation_durations(1:length(list_saccade_aligneddata.trial));
%     stretched_time_to_fixation = stretched_time_to_fixation(1:length(sequence_saccade_aligneddata.trial));
%
%     if strcmpi(a{files(f)}(1:2),'PW')
%         predict_rt = 156;
%     else %for Tobii
%         predict_rt = 138;
%     end
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%-Remove Trials with NaNs--%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %in case it was not implemented before (i.e. during first import)
%
%
%     %---Saccades List trials---%
%     %remove trials with NaNs
%     remove_trials = [];
%     for t = 1:length(list_saccade_aligneddata.trial)
%         if sum(sum(isnan(list_saccade_aligneddata.trial{t}))) >1
%             remove_trials =[remove_trials t];
%         end
%     end
%     list_saccade_aligneddata.time(remove_trials) = [];
%     list_saccade_aligneddata.trial(remove_trials) = [];
%     list_saccade_aligneddata.sampleinfo(remove_trials,:) = [];
%     stretched_fixation_durations(remove_trials) = [];
%
%      %---Saccades Sequence trials---%
%     %remove trials with NaNs
%     remove_trials = [];
%     for t = 1:length(sequence_saccade_aligneddata.trial)
%         if sum(sum(isnan(sequence_saccade_aligneddata.trial{t}))) >1
%             remove_trials =[remove_trials t];
%         end
%     end
%     sequence_saccade_aligneddata.time(remove_trials) = [];
%     sequence_saccade_aligneddata.trial(remove_trials) = [];
%     sequence_saccade_aligneddata.sampleinfo(remove_trials,:) = [];
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%---Calcaulate Baseline (ITI) Power Spectra---%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     %---For ListSQ---%
%     list_low_freq_ITI = lfp_powerspectrum(listsq_ITI_aligneddata,'low','all',time_window);
%     list_high_freq_ITI = lfp_powerspectrum(listsq_ITI_aligneddata,'high','all',time_window);
%     list_low_freq_ITI_BaseLinePower = Seths_BaselinePower(list_low_freq_ITI,baseline_time_interval);
%     list_high_freq_ITI_BaseLinePower = Seths_BaselinePower(list_high_freq_ITI,baseline_time_interval);
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%---Calculate LFP Waveforms, PowerSpectra, and InterSaccadic-Coherence---%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     cfg = [];
%     cfg.channel   = 'all';
%     cfg.trials    = 'all';
%     cfg.padding   = 1;
%
%     %---For List Trials---%
%     %caluculate average LFP waveform for various duration saccades in List task
%     for fdt = 1:size(fixdurthresh,1)
%         fixind = find(stretched_fixation_durations >= fixdurthresh(fdt,1) & stretched_fixation_durations <= fixdurthresh(fdt,2));
%         if length(fixind) > 10 %need at least 10
%             temp_aligned = list_saccade_aligneddata;
%             temp_aligned.trial = temp_aligned.trial(fixind);
%             temp_aligned.time = temp_aligned.time(fixind);
%             temp_aligned.sampleinfo = temp_aligned.sampleinfo(fixind,:);
%             temp_aligned.cfg.trl = temp_aligned.cfg.trl(fixind);
%
%             %Get Average Waveform
%             timelock = ft_timelockanalysis(cfg,temp_aligned);
%             list_average_LFP{fdt,f} = timelock.avg;
%
%             %Get Low and High Frequency Power Spectra
%             freq = lfp_powerspectrum(temp_aligned,'low','all',time_window);
%             %normalize to baseline power by frequency and channel
%             list_low_freq{fdt,f} = Apply_Seths_BaselinePower(freq,list_low_freq_ITI_BaseLinePower,'relchange');
%
%             freq = lfp_powerspectrum(temp_aligned,'high','all',time_window);
%             %normalize to baseline power by frequency and channel
%             list_high_freq{fdt,f} = Apply_Seths_BaselinePower(freq,list_high_freq_ITI_BaseLinePower,'relchange');
%
%             %Get Low and High Frequency Coherence
%             list_low_freq_coh{fdt,f} = lfp_phasecoherence(temp_aligned,'low','all',time_window);
%             list_high_freq_coh{fdt,f} = lfp_phasecoherence(temp_aligned,'high','all',time_window);
%         end
%     end
%
%     cfg = [];
%     cfg.channel   = 'all';
%     cfg.trials    = 'all';
%     cfg.padding   = 1;
%
%     %---For Sequence Trials---%
%     %Remove predictive saccades only want reactive saccades incase there's an
%     %influence of visual stimulus on LFPs or corrective saccades
%     fixind = find(stretched_time_to_fixation < predict_rt);
%     sequence_saccade_aligneddata.trial(fixind) = [];
%     sequence_saccade_aligneddata.time(fixind) = [];
%
%     %Get Average Waveform
%     timelock = ft_timelockanalysis(cfg,sequence_saccade_aligneddata);
%     sequence_average_LFP{f} = timelock.avg;
%
%
%     %Get Low and High Frequency Power Spectra
%     freq = lfp_powerspectrum(sequence_saccade_aligneddata,'low','all',time_window);
%     %normalize to baseline power by frequency and channel
%     sequence_low_freq{f} = Apply_Seths_BaselinePower(freq,list_low_freq_ITI_BaseLinePower,'relchange');
%
%     freq = lfp_powerspectrum(sequence_saccade_aligneddata,'high','all',time_window);
%     %normalize to baseline power by frequency and channel
%     sequence_high_freq{f} = Apply_Seths_BaselinePower(freq,list_high_freq_ITI_BaseLinePower,'relchange');
%
%     %Get Low and High Frequency Coherence
%     sequence_low_freq_coh{f} = lfp_phasecoherence(sequence_saccade_aligneddata,'low','all',time_window);
%     sequence_high_freq_coh{f} = lfp_phasecoherence(sequence_saccade_aligneddata,'high','all',time_window);
% end
% save([data_dir 'ListSQ_saccade_aligned_averaged_across_sessions.mat'],...
%     'list_average_LFP','list_low_freq','list_high_freq','list_low_freq_coh',...
%     'list_high_freq_coh','sequence_average_LFP','sequence_low_freq',...
%     'sequence_high_freq','sequence_low_freq_coh','sequence_high_freq_coh');
%
%
% clear('list_saccade_aligneddata','stretched_fixation_durations','temp_aligned',...
%     'sequence_saccade_aligneddata','stretched_time_to_fixation','listsq_ITI_aligneddata')
%
% %%
% load(['C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\'...
% 'ListSQ_saccade_aligned_averaged_across_sessions.mat']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Combined Data Across Sessions---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%put in structure so we can do statistical analysis later. Not an efficeint
%way by any means to organize data for averaging but it is the most unviersal

%---List Variables---%
all_session_average_waveform = cell(1,length(fixdurthresh));
all_session_low_freq = cell(length(fixdurthresh),length(list_low_freq{1}.time),length(list_low_freq{1}.freq));
all_session_high_freq = cell(length(fixdurthresh),length(list_low_freq{1}.time),length(list_high_freq{1}.freq));
all_session_low_freq_coh = cell(length(fixdurthresh),length(list_low_freq{1}.time),length(list_low_freq{1}.freq));
all_session_high_freq_coh = cell(length(fixdurthresh),length(list_low_freq{1}.time),length(list_high_freq{1}.freq));

%---Sequence Variables---%
all_sequence_averages = [];
sequence_all_session_low_freq = cell(length(sequence_low_freq{1}.time),length(sequence_low_freq{1}.freq));
sequence_all_session_high_freq = cell(length(sequence_high_freq{1}.time),length(sequence_high_freq{1}.freq));
sequence_all_session_low_freq_coh = cell(length(sequence_low_freq{1}.time),length(sequence_low_freq{1}.freq));
sequence_all_session_high_freq_coh = cell(length(sequence_high_freq{1}.time),length(sequence_high_freq{1}.freq));


for f = 1:length(list_average_LFP)
    %---For List---%
    for fdt = 1:length(fixdurthresh)
        all_session_average_waveform{fdt} = [all_session_average_waveform{fdt}; list_average_LFP{fdt,f}];%all "waveforms"
        for chan = 1:size(list_low_freq{fdt,f}.powspctrm,1)
            for t = 1:length(list_high_freq{1}.time)

                %for low frequency
                for freq = 1:length(list_low_freq{1}.freq)
                    all_session_low_freq{fdt,t,freq} =[all_session_low_freq{fdt,t,freq} list_low_freq{fdt,f}.powspctrm(chan,freq,t)];
                    all_session_low_freq_coh{fdt,t,freq} =[all_session_low_freq_coh{fdt,t,freq} list_low_freq_coh{fdt,f}.cohspctrm(chan,freq,t)];
                end


                %for high frequency
                for freq = 1:length(list_high_freq{1}.freq)
                    all_session_high_freq{fdt,t,freq} =[all_session_high_freq{fdt,t,freq} list_high_freq{fdt,f}.powspctrm(chan,freq,t)];
                    all_session_high_freq_coh{fdt,t,freq} =[all_session_high_freq_coh{fdt,t,freq} list_high_freq_coh{fdt,f}.cohspctrm(chan,freq,t)];
                end
            end
        end
    end

    %---For Sequence---%
    all_sequence_averages = [all_sequence_averages; sequence_average_LFP{f}];
    for chan = 1:size(sequence_low_freq{f}.powspctrm,1)
        for t = 1:length(sequence_high_freq{1}.time)

            %for low frequency
            for freq = 1:length(sequence_low_freq{1}.freq)
                sequence_all_session_low_freq{t,freq} =[sequence_all_session_low_freq{t,freq} sequence_low_freq{f}.powspctrm(chan,freq,t)];
                sequence_all_session_low_freq_coh{t,freq} =[sequence_all_session_low_freq_coh{t,freq} sequence_low_freq_coh{f}.cohspctrm(chan,freq,t)];
            end


            %for high frequency
            for freq = 1:length(sequence_high_freq{1}.freq)
                sequence_all_session_high_freq{t,freq} =[sequence_all_session_high_freq{t,freq} sequence_high_freq{f}.powspctrm(chan,freq,t)];
                sequence_all_session_high_freq_coh{t,freq} =[sequence_all_session_high_freq_coh{t,freq} sequence_high_freq_coh{f}.cohspctrm(chan,freq,t)];
            end
        end
    end
end
%%
%---Calculate averages---%
all_session_average_low_freq = cell(1,length(fixdurthresh));
all_session_average_high_freq = cell(1,length(fixdurthresh));
all_session_average_low_freq_coh = cell(1,length(fixdurthresh));
all_session_average_high_freq_coh = cell(1,length(fixdurthresh));

%for list
for fdt = 1:length(fixdurthresh)
    for t = 1:length(list_high_freq{1}.time)

        %for low frequency
        for freq = 1:length(list_low_freq{1}.freq)
            all_session_average_low_freq{fdt}(freq,t) = nanmean(all_session_low_freq{fdt,t,freq});
            all_session_average_low_freq_coh{fdt}(freq,t) = nanmean(all_session_low_freq_coh{fdt,t,freq});
        end

        %for high frequency
        for freq = 1:length(list_high_freq{1}.freq)
            all_session_average_high_freq{fdt}(freq,t) = nanmean(all_session_high_freq{fdt,t,freq});
            all_session_average_high_freq_coh{fdt}(freq,t) = nanmean(all_session_high_freq_coh{fdt,t,freq});
        end
    end
end

%for sequence
sequence_all_session_average_low_freq = NaN(length(sequence_low_freq{1}.freq),length(sequence_low_freq{1}.time));
sequence_all_session_average_high_freq = NaN(length(sequence_high_freq{1}.freq),length(sequence_low_freq{1}.time));
sequence_all_session_average_low_freq_coh = NaN(length(sequence_low_freq{1}.freq),length(sequence_low_freq{1}.time));
sequence_all_session_average_high_freq_coh = NaN(length(sequence_high_freq{1}.freq),length(sequence_low_freq{1}.time));
for t = 1:length(sequence_low_freq{1}.time)

    %for low frequency
    for freq = 1:length(sequence_low_freq{1}.freq)
        sequence_all_session_average_low_freq(freq,t) = nanmean(sequence_all_session_low_freq{t,freq});
        sequence_all_session_average_low_freq_coh(freq,t) = nanmean(sequence_all_session_low_freq_coh{t,freq});
    end

    %for high frequency
    for freq = 1:length(sequence_high_freq{1}.freq)
        sequence_all_session_average_high_freq(freq,t) = nanmean(sequence_all_session_high_freq{t,freq});
        sequence_all_session_average_high_freq_coh(freq,t) = nanmean(sequence_all_session_high_freq_coh{t,freq});
    end
end

%%
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Summary Figures\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot average LFP waveforms---%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = (-2048:2048)/1000;
smval = [];
colors ={'cyan','green','blue','magenta','black'};

all_session_averages = cell(1,size(fixdurthresh,1));
all_sequence_averages = [];
for f = 1:size(list_average_LFP,2)
    for fdt = 1:size(fixdurthresh,1)
        all_session_averages{fdt} = [all_session_averages{fdt}; list_average_LFP{fdt,f}];
    end
    all_sequence_averages = [all_sequence_averages; sequence_average_LFP{f}];
end

figure
hold on
lg = cell(1,size(fixdurthresh,1));
for fdt = 1:size(fixdurthresh,1)
    dofill(time,all_session_averages{fdt},colors{fdt},smval);
    lg{fdt} = ['Fix durs: ' num2str(fixdurthresh(fdt,1)) '-' num2str( fixdurthresh(fdt,2))];
end
dofill(time,all_sequence_averages,'red',smval);
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.5 0.8])
xlabel('Time From Saccade Start')
ylabel('LFP (uV)')
legend([lg {'Reactive Sequence'}])
title(sprintf(['Saccade Locked LFPs across tasks and Fixation durations n_{electrode} = '...
    num2str(size(all_sequence_averages,1))]));
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
grid(gca,'minor')

save_and_close_fig(figure_dir,'List_by_fixdur_and_Sequence-Waveforms')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Low Frequency Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_time = list_low_freq{1}.time;
low_freq = list_low_freq{1}.freq;
high_freq = list_high_freq{1}.freq;


yl = NaN(6,2);

figure

%---List Saccades---%
for fdt = 1:length(fixdurthresh)
    subplot(2,3,fdt)
    imagesc(freq_time,low_freq,all_session_average_low_freq{fdt})
    ylabel('Frequency (Hz)')
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(sprintf(['List Fixdur = ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2)) ' ms']))
end

%---Sequene Saccades---%
subplot(2,3,6)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq)
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(6,:) = caxis;
title('Sequence')

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(low_freq) max(low_freq)],'w--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

subtitle('List and Sequence Saccades: Low Frequencey Power Spectral Analysis')
save_and_close_fig(figure_dir,'List_by_fixdur_and_Sequence-Low Frequency Power Spectrum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot high Frequency Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_time = list_high_freq{1}.time;
high_freq = list_high_freq{1}.freq;
high_freq = list_high_freq{1}.freq;


yl = NaN(6,2);

figure

%---List Saccades---%
for fdt = 1:length(fixdurthresh)
    subplot(2,3,fdt)
    imagesc(freq_time,high_freq,all_session_average_high_freq{fdt})
    ylabel('Frequency (Hz)')
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(sprintf(['List Fixdur = ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2)) ' ms']))
end

%---Sequene Saccades---%
subplot(2,3,6)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq)
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(6,:) = caxis;
title('Sequence')

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(high_freq) max(high_freq)],'w--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

subtitle('List and Sequence Saccades: high Frequencey Power Spectral Analysis')
save_and_close_fig(figure_dir,'List_by_fixdur_and_Sequence-High Frequency Power Spectrum')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Low Frequency Coherence Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yl = NaN(6,2);

figure

%---List Saccades---%
for fdt = 1:length(fixdurthresh)
    subplot(2,3,fdt)
    imagesc(freq_time,low_freq,all_session_average_low_freq_coh{fdt})
    ylabel('Frequency (Hz)')
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(sprintf(['List Fixdur = ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2)) ' ms']))
end

%---Sequene Saccades---%
subplot(2,3,6)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq_coh)
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(6,:) = caxis;
title('Sequence')

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(high_freq) max(high_freq)],'w--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

subtitle('List and Sequence Saccades: Low Frequencey Coherence Spectral Analysis')
save_and_close_fig(figure_dir,'List_by_fixdur_and_Sequence-Low Frequency Coherence Spectrum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot high Frequency Coherence Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yl = NaN(6,2);

figure

%---List Saccades---%
for fdt = 1:length(fixdurthresh)
    subplot(2,3,fdt)
    imagesc(freq_time,high_freq,all_session_average_high_freq_coh{fdt})
    ylabel('Frequency (Hz)')
    xlim([-0.75 0.75])
    xlabel('Time from Saccade (sec)')
    yl(fdt,:) = caxis;
    title(sprintf(['List Fixdur = ' num2str(fixdurthresh(fdt,1)) '-' num2str(fixdurthresh(fdt,2)) ' ms']))
end

%---Sequene Saccades---%
subplot(2,3,6)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq_coh)
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(6,:) = caxis;
title('Sequence')

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:6
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(high_freq) max(high_freq)],'w--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

subtitle('List and Sequence Saccades: high Frequencey Coherence Spectral Analysis')
save_and_close_fig(figure_dir,'List_by_fixdur_and_Sequence-High Frequency Coherence Spectrum')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Frequency Statistical Anlaysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = .150;
pre_time = find(freq_time >= -twin & freq_time  < 0);
post_time = find(freq_time >= 0 & freq_time < twin);
band  = [4 12; %theta
    12 16; %alpha
    16 30; %beta
    30 50; %low gamma
    50 120]; %high gamma

freq_pre_post_ps = NaN(1,size(band,1));
freq_means = NaN(2,size(band,1));
freq_stds = NaN(2,size(band,1));


seq_freq_pre_post_ps = NaN(1,size(band,1));
seq_freq_means = NaN(2,size(band,1));
seq_freq_stds = NaN(2,size(band,1));

%for frequency spectrum analysis
figure
for b = 1:size(band,1);
    
    %---For List---%
    if b <= 3 %for the low frequencies
        freq_ind = find(low_freq >= band(b,1) & low_freq <= band(b,2));
        
        subplot(2,3,b)
        hold on
        for fdt = 1:length(fixdurthresh)
            %grab the data
            temp = cell(t,length(freq_ind));
            for f = 1:length(freq_ind)
                for t = 1:size(all_session_low_freq,2)
                    temp{t,f} = all_session_low_freq{fdt,t,freq_ind(f)};
                end
            end
            %average across the frequencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(all_session_low_freq,2));
            for t = 1:size(all_session_low_freq,2)
               temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
            end
            dofill(freq_time,temp2,colors{fdt},smval);            
        end
        
        %do for last fdt so 500-750 ms since closet to sequence
        pre_window = nanmean(temp2(:,pre_time),2);
        freq_means(1,b) = nanmean(pre_window);
        freq_stds(1,b) = nanstd(pre_window);
        
        post_window = nanmean(temp2(:,post_time));
        freq_means(2,b) = nanmean(post_window);
        freq_stds(2,b) = nanstd(post_window);
        
        [~,freq_pre_post_ps] = ttest2(post_window,pre_window);

        
    else %for the high frequencies
        freq_ind = find(high_freq >= band(b,1) & high_freq <= band(b,2));

        subplot(2,3,b)
        hold on
        for fdt = 1:length(fixdurthresh)
            %grab the data
            temp = cell(t,length(freq_ind));
            for f = 1:length(freq_ind)
                for t = 1:size(all_session_high_freq,2)
                    temp{t,f} = all_session_high_freq{fdt,t,freq_ind(f)};
                end
            end
            %average across the frequencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(all_session_high_freq,2));
            for t = 1:size(all_session_high_freq,2)
               temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
            end
            dofill(freq_time,temp2,colors{fdt},smval); 
            
            pre_window = nanmean(temp2(:,pre_time),2);
            freq_means(1,b) = nanmean(pre_window);
            freq_stds(1,b) = nanstd(pre_window);
            
            post_window = nanmean(temp2(:,post_time));
            freq_means(2,b) = nanmean(post_window);
            freq_stds(2,b) = nanstd(post_window);
            
            [~,freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
        end
    end
    
     %---For Sequence---%
    if b <= 3 %for the low frequencies
        freq_ind = find(low_freq >= band(b,1) & low_freq <= band(b,2));
        
        temp = sequence_all_session_low_freq(:,freq_ind);
        
        subplot(2,3,b)
        hold on
        %average across the frequencies in the band but keep
        %electrode-session seperate
        temp2 = NaN(size(all_session_low_freq{1},2),size(temp,1));
        for t = 1:size(temp2,1)
            temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
        end
        dofill(freq_time,temp2,'red',smval);
        
    else %for the high frequencies
          freq_ind = find(high_freq >= band(b,1) & high_freq <= band(b,2));
        
        temp = sequence_all_session_high_freq(:,freq_ind);
        
        subplot(2,3,b)
        hold on
        %average across the frequencies in the band but keep
        %electrode-session seperate
        temp2 = NaN(size(all_session_high_freq{1},2),size(temp,1));
        for t = 1:size(temp2,1)
            temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
        end
        dofill(freq_time,temp2,'red',smval);
    end
      pre_window = nanmean(temp2(:,pre_time),2);
        seq_freq_means(1,b) = nanmean(pre_window);
        seq_freq_stds(1,b) = nanstd(pre_window);
        
        post_window = nanmean(temp2(:,post_time));
        seq_freq_means(2,b) = nanmean(post_window);
        seq_freq_stds(2,b) = nanstd(post_window);
        
        [~,freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
    
    
end

band_str = {'Theta','Alpha','Beta','Low Gamma','High Gamma'};
nsamp = sqrt(size(temp2,1));
for sb = 1:5
    subplot(2,3,sb)
    yl = ylim;
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-0.75 0.75])
    xlabel('Time From Saccade Start')
    ylabel('Relative Power')
    title([band_str{sb} ': ' num2str(band(sb,1)) '-'  num2str(band(sb,2)) ' Hz'])
    
    if sb == 5
       legend([lg,{'Sequence'}]) 
    end
end

%bar graph for power before and after saccade
subplot(4,3,9) %for list only 500-600 ms
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[freq_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[freq_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[freq_means(1,b) freq_means(2,b)],[freq_stds(1,b) freq_stds(2,b)]./nsamp)
    if freq_pre_post_ps(b) < 0.05
       plot(b,0.2,'*k') 
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Relative Power')
title('List Fix durs 500-700')

subplot(4,3,12) %for squence
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[seq_freq_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[seq_freq_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[seq_freq_means(1,b) seq_freq_means(2,b)],[seq_freq_stds(1,b) seq_freq_stds(2,b)]./nsamp)
    if seq_freq_pre_post_ps(b) < 0.25
       plot(b,0.5,'*k') 
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Relative Power')
title('Sequence')

save_and_close_fig(figure_dir,'List_by_fixdur_and_Sequence-Line Plot Frequency Power Spectrum')
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Coherence Statistical Anlaysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


coh_freq_pre_post_ps = NaN(1,size(band,1));
coh_freq_means = NaN(2,size(band,1));
coh_freq_stds = NaN(2,size(band,1));


seq_coh_freq_pre_post_ps = NaN(1,size(band,1));
seq_coh_freq_means = NaN(2,size(band,1));
seq_coh_freq_stds = NaN(2,size(band,1));

%for frequency spectrum analysis
figure
for b = 1:size(band,1);
    
    %---For List---%
    if b <= 3 %for the low frequencies
        freq_ind = find(low_freq >= band(b,1) & low_freq <= band(b,2));
        
        subplot(2,3,b)
        hold on
        for fdt = 1:length(fixdurthresh)
            %grab the data
            temp = cell(t,length(freq_ind));
            for f = 1:length(freq_ind)
                for t = 1:size(all_session_low_freq,2)
                    temp{t,f} = all_session_low_freq_coh{fdt,t,freq_ind(f)};
                end
            end
            %average across the frequencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(all_session_low_freq,2));
            for t = 1:size(all_session_low_freq,2)
                temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
            end
            dofill(freq_time,temp2-min(nanmean(temp2)),colors{fdt},smval);
        end
        
        %do for last fdt so 500-750 ms since closet to sequence
        pre_window = nanmean(temp2(:,pre_time),2);
        coh_freq_means(1,b) = nanmean(pre_window);
        coh_freq_stds(1,b) = nanstd(pre_window);
        
        post_window = nanmean(temp2(:,post_time));
        coh_freq_means(2,b) = nanmean(post_window);
        coh_freq_stds(2,b) = nanstd(post_window);
        
        [~,coh_freq_pre_post_ps] = ttest2(post_window,pre_window);
        
        
    else %for the high frequencies
        freq_ind = find(high_freq >= band(b,1) & high_freq <= band(b,2));
        
        subplot(2,3,b)
        hold on
        for fdt = 1:length(fixdurthresh)
            %grab the data
            temp = cell(t,length(freq_ind));
            for f = 1:length(freq_ind)
                for t = 1:size(all_session_high_freq,2)
                    temp{t,f} = all_session_high_freq_coh{fdt,t,freq_ind(f)};
                end
            end
            %average across the frequencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(all_session_high_freq,2));
            for t = 1:size(all_session_high_freq,2)
                temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
            end
            dofill(freq_time,temp2-min(nanmean(temp2)),colors{fdt},smval);
            
            pre_window = nanmean(temp2(:,pre_time),2);
            coh_freq_means(1,b) = nanmean(pre_window);
            coh_freq_stds(1,b) = nanstd(pre_window);
            
            post_window = nanmean(temp2(:,post_time));
            coh_freq_means(2,b) = nanmean(post_window);
            coh_freq_stds(2,b) = nanstd(post_window);
            
            [~,coh_freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
        end
    end
    
    %---For Sequence---%
    if b <= 3 %for the low frequencies
        freq_ind = find(low_freq >= band(b,1) & low_freq <= band(b,2));
        
        temp = sequence_all_session_low_freq_coh(:,freq_ind);
        
        subplot(2,3,b)
        hold on
        %average across the frequencies in the band but keep
        %electrode-session seperate
        temp2 = NaN(size(all_session_low_freq_coh{1},2),size(temp,1));
        for t = 1:size(temp2,1)
            temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
        end
        dofill(freq_time,temp2-min(nanmean(temp2)),'red',smval);
        
    else %for the high frequencies
        freq_ind = find(high_freq >= band(b,1) & high_freq <= band(b,2));
        
        temp = sequence_all_session_high_freq_coh(:,freq_ind);
        
        subplot(2,3,b)
        hold on
        %average across the frequencies in the band but keep
        %electrode-session seperate
        temp2 = NaN(size(all_session_high_freq_coh{1},2),size(temp,1));
        for t = 1:size(temp2,1)
            temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
        end
        dofill(freq_time,temp2-min(nanmean(temp2)),'red',smval);
    end
    pre_window = nanmean(temp2(:,pre_time),2);
    coh_seq_freq_means(1,b) = nanmean(pre_window);
    coh_seq_freq_stds(1,b) = nanstd(pre_window);
    
    post_window = nanmean(temp2(:,post_time));
    coh_seq_freq_means(2,b) = nanmean(post_window);
    coh_seq_freq_stds(2,b) = nanstd(post_window);
    
    [~,coh_freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
    
    
end

band_str = {'Theta','Alpha','Beta','Low Gamma','High Gamma'};
nsamp = sqrt(size(temp2,1));
for sb = 1:5
    subplot(2,3,sb)
    yl = ylim;
    if yl(1) < 0
        ylim([0 yl(2)])
        yl(1) = 0;
    end
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-0.75 0.75])
    xlabel('Time From Saccade Start')
    ylabel('Baseline-Corrected Coherence')
    title([band_str{sb} ': ' num2str(band(sb,1)) '-'  num2str(band(sb,2)) ' Hz'])
    
    if sb == 5
        legend([lg,{'Sequence'}])
        xlim([-0.2 .2])
    end
end

%bar graph for power before and after saccade
subplot(4,3,9) %for list only 500-600 ms
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[coh_freq_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[coh_freq_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[coh_freq_means(1,b) coh_freq_means(2,b)],[coh_freq_stds(1,b) coh_freq_stds(2,b)]./nsamp)
    if freq_pre_post_ps(b) < 0.05
        plot(b,0.2,'*k')
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Coherence')
title('List Fix durs 500-700')

subplot(4,3,12) %for squence
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[coh_seq_freq_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[coh_seq_freq_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[coh_seq_freq_means(1,b) coh_seq_freq_means(2,b)],[coh_seq_freq_stds(1,b) coh_seq_freq_stds(2,b)]./nsamp)
    if seq_freq_pre_post_ps(b) < 0.25
        plot(b,0.5,'*k')
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Coherence')
title('Sequence')
save_and_close_fig(figure_dir,'List_by_fixdur_and_Sequence-Line Plot Coherence Power Spectrum')

