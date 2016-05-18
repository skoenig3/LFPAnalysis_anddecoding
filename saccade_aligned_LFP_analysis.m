function saccade_aligned_LFP_analysis(session_data,data_dir,figure_dir)
%written by Seth Konig 4/28/16
%code analyses power spectral properties of LFPs aligned to eye movements
%across List, Sequecne, and Covert attention sessions.

task_file = get_task_data(session_data,'ListSQ');
if isempty(task_file)
    return
end

if exist([data_dir task_file(1:8) '-preproccessed_saccade_aligned_LFPs.mat'],'file')
    load([data_dir task_file(1:8) '-preproccessed_saccade_aligned_LFPs.mat'])
else
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---remove any trials with NaNs, code can't process them well--%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for list saccades
remove_trials = [];
for t = 1:length(list_saccade_aligneddata.trial)
    if sum(sum(isnan(list_saccade_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
list_saccade_aligneddata.time(remove_trials) = [];
list_saccade_aligneddata.trial(remove_trials) = [];
list_saccade_aligneddata.sampleinfo(remove_trials,:) = [];

%for sequence saccades
remove_trials = [];
for t = 1:length(sequence_saccade_aligneddata.trial)
    if sum(sum(isnan(sequence_saccade_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
sequence_saccade_aligneddata.time(remove_trials) = [];
sequence_saccade_aligneddata.trial(remove_trials) = [];
sequence_saccade_aligneddata.sampleinfo(remove_trials,:) = [];

%for listsq ITI period
remove_trials = [];
for t = 1:length(listsq_ITI_aligneddata.trial)
    if sum(sum(isnan(listsq_ITI_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
listsq_ITI_aligneddata.time(remove_trials) = [];
listsq_ITI_aligneddata.trial(remove_trials) = [];
listsq_ITI_aligneddata.sampleinfo(remove_trials,:) = [];

%for cvtnew ITI period
remove_trials = [];
for t = 1:length(cvtnew_ITI_aligneddata.trial)
    if sum(sum(isnan(cvtnew_ITI_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
cvtnew_ITI_aligneddata.time(remove_trials) = [];
cvtnew_ITI_aligneddata.trial(remove_trials) = [];
cvtnew_ITI_aligneddata.sampleinfo(remove_trials,:) = [];

%for cvtnew fixation
remove_trials = [];
for t = 1:length(cvtfix_aligneddata.trial)
    if sum(sum(isnan(cvtfix_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
cvtfix_aligneddata.time(remove_trials) = [];
cvtfix_aligneddata.trial(remove_trials) = [];
cvtfix_aligneddata.sampleinfo(remove_trials,:) = [];

%for cvtnew dot aligned data
remove_trials = [];
for t = 1:length(cvtfix_aligneddata.trial)
    if sum(sum(isnan(dot_aligneddata.trial{t}))) >1
        remove_trials =[remove_trials t];
    end
end
dot_aligneddata.time(remove_trials) = [];
dot_aligneddata.trial(remove_trials) = [];
dot_aligneddata.sampleinfo(remove_trials,:) = [];

baseline_time_interval = [0.25 0.75]; %in ITI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Calculate Baseline Power---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfgfrq = [];
cfgfrq.baseline = 'no';
cfgfrq.baselinetype = [];
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';
time_window = -0.75:.01:1.5;

%---For ListSQ---%
list_low_freq_ITI = lfp_powerspectrum(listsq_ITI_aligneddata,'low','all',time_window);
list_high_freq_ITI = lfp_powerspectrum(listsq_ITI_aligneddata,'high','all',time_window);
list_low_freq_ITI_BaseLinePower = Seths_BaselinePower(list_low_freq_ITI,baseline_time_interval);
list_high_freq_ITI_BaseLinePower = Seths_BaselinePower(list_high_freq_ITI,baseline_time_interval);

%---For CVTNEW---%
cvtnew_low_freq_ITI = lfp_powerspectrum(cvtnew_ITI_aligneddata,'low','all',time_window);
cvtnew_high_freq_ITI = lfp_powerspectrum(cvtnew_ITI_aligneddata,'high','all',time_window);
cvtnew_low_freq_ITI_BaseLinePower = Seths_BaselinePower(cvtnew_low_freq_ITI,baseline_time_interval);
cvtnew_high_freq_ITI_BaseLinePower = Seths_BaselinePower(cvtnew_high_freq_ITI,baseline_time_interval);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%---Plot All LFP "Waveforms"---%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Processing Average LFP "Wavforms" for ' task_file(1:8)])

average_LFP_waveforms = cell(1,4); %average waveforms
% list, sequence, cvtnew fixonset, cvtnew dot on

cfg = [];
cfg.channel   = 'all';
cfg.trials    = 'all';
cfg.padding   = 1;

yl = NaN(4,2);
figure

%---Saccades List trials---%
timelock = ft_timelockanalysis(cfg,list_saccade_aligneddata);
average_LFP_waveforms{1} = timelock;
subplot(2,2,1)
plot(timelock.time,timelock.avg')
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(1,:) = ylim;
title(sprintf(['List Trials n_{saccades} = ' num2str(size(list_saccade_aligneddata.sampleinfo,1))]))

%---Saccades Sequence trials---%
timelock = ft_timelockanalysis(cfg,sequence_saccade_aligneddata);
average_LFP_waveforms{2} = timelock;
subplot(2,2,2)
plot(timelock.time,timelock.avg')
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(2,:) = ylim;
title(sprintf(['Sequence Trials n_{saccades} = ' num2str(size(sequence_saccade_aligneddata.sampleinfo,1))]))

%---CVTNEW Fixation aligned trials---%
timelock = ft_timelockanalysis(cfg,cvtfix_aligneddata);
average_LFP_waveforms{3} = timelock;
subplot(2,2,3)
plot(timelock.time,timelock.avg')
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Fixation (sec)')
ylabel('Voltage (uV)')
yl(3,:) = ylim;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---CVTNEW dot aligned trials---%
timelock = ft_timelockanalysis(cfg,dot_aligneddata);
average_LFP_waveforms{4} = timelock;
subplot(2,2,4)
plot(timelock.time,timelock.avg')
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Dot on (sec)')
ylabel('Voltage (uV)')
yl(4,:) = ylim;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:4
    subplot(2,2,sb)
    hold on
    line([0 0],[ymin ymax],'Color','k','LineStyle','--')
    hold off
    set(gca,'XMinorTick','on','YMinorTick','on')
    grid on
    grid(gca,'minor')
    ylim([ymin ymax])
end
save_and_close_fig(figure_dir,[task_file(1:8) '-Across_tasks_saccade_aligned_waveforms'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%---Plot All Low Freq Power Spectrum---%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Processing Low Frequency Power Spectrum for ' task_file(1:8)])

cfgfrq = [];
cfgfrq.baseline = [-0.75 -0.5];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';
cfgfrq.colormap     = colormap('jet');

low_freq_data = cell(1,4);
yl = NaN(4,2);

figure

%---Saccades List trials---%
freq = lfp_powerspectrum(list_saccade_aligneddata,'low','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,list_low_freq_ITI_BaseLinePower,'relchange');
low_freq_data{1} = freq;
subplot(2,2,1)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(1,:) = caxis;
title(sprintf(['List Trials n_{saccades} = ' num2str(size(list_saccade_aligneddata.sampleinfo,1))]))

%---Saccades Sequence Trials---%
freq = lfp_powerspectrum(sequence_saccade_aligneddata,'low','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,list_low_freq_ITI_BaseLinePower,'relchange');
low_freq_data{2} = freq;
subplot(2,2,2)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(2,:) = caxis;
title(sprintf(['List Trials n_{saccades} = ' num2str(size(sequence_saccade_aligneddata.sampleinfo,1))]))

%---CVTNEW Fixation aligned trials---%
freq = lfp_powerspectrum(cvtfix_aligneddata,'low','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,cvtnew_low_freq_ITI_BaseLinePower,'relchange');
low_freq_data{3} = freq;
subplot(2,2,3)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Fixation (sec)')
ylabel('Voltage (uV)')
yl(3,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---CVTNEW dot aligned trials---%
freq = lfp_powerspectrum(dot_aligneddata,'low','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,cvtnew_low_freq_ITI_BaseLinePower,'relchange');
low_freq_data{4} = freq;
subplot(2,2,4)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(4,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:2
    subplot(2,2,sb)
    hold on
    plot([0 0],[min(freq.freq) max(freq.freq)],'w--')
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,[task_file(1:8) '-Across_tasks_saccade_aligned_low_freq_powerspectrum'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---Plot All high Freq Power Spectrum---%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Processing High Frequency Power Spectrum for ' task_file(1:8)])

cfgfrq = [];
cfgfrq.baseline = [-0.75 -0.5];
cfgfrq.baselinetype = 'absolute';
cfgfrq.maskstyle    = 'saturation';
cfgfrq.zparam       = 'powspctrm';
cfgfrq.colormap     = colormap('jet');

high_freq_data = cell(1,4);
yl = NaN(4,2);

figure

%---Saccades List trials---%
freq = lfp_powerspectrum(list_saccade_aligneddata,'high','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,list_high_freq_ITI_BaseLinePower,'relchange');
high_freq_data{1} = freq;
subplot(2,2,1)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(1,:) = caxis;
title(sprintf(['List Trials n_{saccades} = ' num2str(size(list_saccade_aligneddata.sampleinfo,1))]))

%---Saccades Sequence Trials---%
freq = lfp_powerspectrum(sequence_saccade_aligneddata,'high','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,list_high_freq_ITI_BaseLinePower,'relchange');
high_freq_data{2} = freq;
subplot(2,2,2)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(2,:) = caxis;
title(sprintf(['List Trials n_{saccades} = ' num2str(size(sequence_saccade_aligneddata.sampleinfo,1))]))

%---CVTNEW Fixation aligned trials---%
freq = lfp_powerspectrum(cvtfix_aligneddata,'high','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,cvtnew_high_freq_ITI_BaseLinePower,'relchange');
high_freq_data{3} = freq;
subplot(2,2,3)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Fixation (sec)')
ylabel('Voltage (uV)')
yl(3,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---CVTNEW dot aligned trials---%
freq = lfp_powerspectrum(dot_aligneddata,'high','all',time_window);
%normalize to baseline power by frequency and channel
freq = Apply_Seths_BaselinePower(freq,cvtnew_high_freq_ITI_BaseLinePower,'relchange');
high_freq_data{4} = freq;
subplot(2,2,4)
ft_singleplotTFR(cfgfrq, freq);
xlim([-0.75 0.75])
set(gca,'Xtick',[-.750:.250:.750])
xlabel('Time from Saccade (sec)')
ylabel('Voltage (uV)')
yl(4,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:4
    subplot(2,2,sb)
    hold on
    plot([0 0],[min(freq.freq) max(freq.freq)],'w--')
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,[task_file(1:8) '-Across_tasks_saccade_aligned_high_freq_powerspectrum'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%---Plot All Low Freq Coherence---%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Processing Low Frequency Intersaccadic Coherence for ' task_file(1:8)])

low_freq_coherence_data = cell(1,4);
yl = NaN(4,2);

figure    


%---Saccades List trials---%
stat = lfp_phasecoherence(list_saccade_aligneddata,'low','all',time_window);
low_freq_coherence_data{1} = stat;
subplot(2,2,1)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title(sprintf(['List Trials n_{saccades} = ' num2str(size(list_saccade_aligneddata.sampleinfo,1))]))

%---Saccades Sequence Trials---%
stat = lfp_phasecoherence(sequence_saccade_aligneddata,'low','all',time_window);
low_freq_coherence_data{2} = stat;
subplot(2,2,2)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title(sprintf(['Sequence Trials n_{saccades} = ' num2str(size(sequence_saccade_aligneddata.sampleinfo,1))]))

%---CVTNEW Fixation aligned trials---%
stat = lfp_phasecoherence(cvtfix_aligneddata,'low','all',time_window);
low_freq_coherence_data{3} = stat;
subplot(2,2,3)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Fixation (sec)')
yl(3,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---CVTNEW dot aligned trials---%
stat = lfp_phasecoherence(dot_aligneddata,'low','all',time_window);
low_freq_coherence_data{4} = stat;
subplot(2,2,4)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Dot on (sec)')
yl(4,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:4
    subplot(2,2,sb)
    hold on
    plot([0 0],[min(stat.freq) max(stat.freq)],'w--')
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,[task_file(1:8) '-Across_tasks_saccade_aligned_low_freq_coherence'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%---Plot All high Freq Coherence---%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

high_freq_coherence_data = cell(1,4);
yl = NaN(4,2);

figure    


%---Saccades List trials---%
stat = lfp_phasecoherence(list_saccade_aligneddata,'high','all',time_window);
high_freq_coherence_data{1} = stat;
subplot(2,2,1)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title(sprintf(['List Trials n_{saccades} = ' num2str(size(list_saccade_aligneddata.sampleinfo,1))]))

%---Saccades Sequence Trials---%
stat = lfp_phasecoherence(sequence_saccade_aligneddata,'high','all',time_window);
high_freq_coherence_data{2} = stat;
subplot(2,2,2)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title(sprintf(['Sequence Trials n_{saccades} = ' num2str(size(sequence_saccade_aligneddata.sampleinfo,1))]))

%---CVTNEW Fixation aligned trials---%
stat = lfp_phasecoherence(cvtfix_aligneddata,'high','all',time_window);
high_freq_coherence_data{3} = stat;
subplot(2,2,3)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Fixation (sec)')
yl(3,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---CVTNEW dot aligned trials---%
stat = lfp_phasecoherence(dot_aligneddata,'high','all',time_window);
high_freq_coherence_data{4} = stat;
subplot(2,2,4)
imagesc(stat.time,stat.freq,reshape(mean(abs(stat.cohspctrm)),size(stat.cohspctrm,2),size(stat.cohspctrm,3)))
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Dot on (sec)')
yl(4,:) = caxis;
title(sprintf(['CVTNEW n_{trials} = ' num2str(size(cvtfix_aligneddata.sampleinfo,1))]))

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:4
    subplot(2,2,sb)
    hold on
    plot([0 0],[min(stat.freq) max(stat.freq)],'w--')
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

save_and_close_fig(figure_dir,[task_file(1:8) '-Across_tasks_saccade_aligned_high_freq_powerspectrum'])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Save the Data---%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
save([data_dir task_file(1:8) '-_saccade_aligned_LFPs_frequency_analysis.mat'],...
    'average_LFP_waveforms','low_freq_data','high_freq_data','low_freq_coherence_data',...
    'high_freq_coherence_data');
end
