% written Seth Konig 4/29/16
% grab analyzed data across all sessions and monkeys and get the averages
% across all sessions and do some list_low_freq_coh{1}istical analysis

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\';

%fixation durations by category minimum already set at 100 ms
Fs = 1000; %sanmpling frequence
smval = []; %no smoothing but want dofill plots
time_window = -0.75:.01:1.5; %what window to consider eye data over
baseline_time_interval = [0.25 0.75]; %in ITI

a = what(data_dir);
a = a.mat;

files = [];
for f = 1:size(a,1)
    if ~isempty(strfind(a{nvr,f},'preproccessed_saccade_aligned_LFPs'))
        files = [f files];
    end
end

%---Preallocate structure for variables---%
list_average_LFP = cell(2,length(files));
list_low_freq = cell(2,length(files));
list_high_freq = cell(2,length(files));
list_low_freq_coh = cell(2,length(files));
list_high_freq_coh = cell(2,length(files));

sequence_average_LFP = cell(1,length(files));
for f = 1:length(files)
    disp([ 'Processing file ' a{files(f)}(1:8)])
    load([data_dir a{files(f)}],'list_saccade_aligneddata','stretched_novel_vs_repeat','listsq_ITI_aligneddata');
    
    %in case it was not implemented before (i.e. during first import)
    stretched_novel_vs_repeat = stretched_novel_vs_repeat(1:length(list_saccade_aligneddata.trial));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%-Remove Trials with NaNs--%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %in case it was not implemented before (i.e. during first import)
    %---Saccades List trials---%
    %remove trials with NaNs
    remove_trials = [];
    for t = 1:length(list_saccade_aligneddata.trial)
        if sum(sum(isnan(list_saccade_aligneddata.trial{t}))) >1
            remove_trials =[remove_trials t];
        end
    end
    list_saccade_aligneddata.time(remove_trials) = [];
    list_saccade_aligneddata.trial(remove_trials) = [];
    list_saccade_aligneddata.sampleinfo(remove_trials,:) = [];
    stretched_novel_vs_repeat(remove_trials) = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calcaulate Baseline (ITI) Power Spectra---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %---For ListSQ---%
    list_low_freq_ITI = lfp_powerspectrum(listsq_ITI_aligneddata,'low','all',time_window);
    list_high_freq_ITI = lfp_powerspectrum(listsq_ITI_aligneddata,'high','all',time_window);
    list_low_freq_ITI_BaseLinePower = Seths_BaselinePower(list_low_freq_ITI,baseline_time_interval);
    list_high_freq_ITI_BaseLinePower = Seths_BaselinePower(list_high_freq_ITI,baseline_time_interval);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%---Calculate LFP Waveforms, PowerSpectra, and InterSaccadic-Coherence---%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cfg = [];
    cfg.channel   = 'all';
    cfg.trials    = 'all';
    cfg.padding   = 1;
    
    %---For List Trials---%
    %caluculate average LFP waveform for various duration saccades in List task
    for nvr = 1:2
        fixind = find(stretched_novel_vs_repeat == nvr);
        temp_aligned = list_saccade_aligneddata;
        temp_aligned.trial = temp_aligned.trial(fixind);
        temp_aligned.time = temp_aligned.time(fixind);
        temp_aligned.sampleinfo = temp_aligned.sampleinfo(fixind,:);
        temp_aligned.cfg.trl = temp_aligned.cfg.trl(fixind);
        
        %Get Average Waveform
        timelock = ft_timelockanalysis(cfg,temp_aligned);
        list_average_LFP{nvr,f} = timelock.avg;
        
        %Get Low and High Frequency Power Spectra
        freq = lfp_powerspectrum(temp_aligned,'low','all',time_window);
        %normalize to baseline power by frequency and channel
        list_low_freq{nvr,f} = Apply_Seths_BaselinePower(freq,list_low_freq_ITI_BaseLinePower,'relchange');
        
        freq = lfp_powerspectrum(temp_aligned,'high','all',time_window);
        %normalize to baseline power by frequency and channel
        list_high_freq{nvr,f} = Apply_Seths_BaselinePower(freq,list_high_freq_ITI_BaseLinePower,'relchange');
        
        %Get Low and High Frequency Coherence
        list_low_freq_coh{nvr,f} = lfp_phasecoherence(temp_aligned,'low','all',time_window);
        list_high_freq_coh{nvr,f} = lfp_phasecoherence(temp_aligned,'high','all',time_window);
    end
end
save([data_dir 'ListSQ_Novel_repeat_saccade_aligned_averaged_across_sessions.mat'],...
    'list_average_LFP','list_low_freq','list_high_freq','list_low_freq_coh',...
    'list_high_freq_coh');


clear('list_saccade_aligneddata','stretched_novel_vs_repeat','temp_aligned','listsq_ITI_aligneddata')
%%
load(['C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\'...
'ListSQ_Novel_repeat_saccade_aligned_averaged_across_sessions.mat']);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Combined Data Across Sessions---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%put in structure so we can do statistical analysis later. Not an efficeint
%way by any means to organize data for averaging but it is the most unviersal
all_session_average_waveform = cell(1,2);
all_session_low_freq = cell(2,length(list_low_freq{1}.time),length(list_low_freq{1}.freq));
all_session_high_freq = cell(2,length(list_low_freq{1}.time),length(list_high_freq{1}.freq));
all_session_low_freq_coh = cell(2,length(list_low_freq{1}.time),length(list_low_freq{1}.freq));
all_session_high_freq_coh = cell(2,length(list_low_freq{1}.time),length(list_high_freq{1}.freq));
for f = 1:length(list_average_LFP)
    for nvr = 1:2
        all_session_average_waveform{nvr} = [all_session_average_waveform{nvr}; list_average_LFP{nvr,f}];%all "waveforms"
        for chan = 1:size(list_low_freq{nvr,f},1)
            for t = 1:length(list_high_freq{1}.time)
                
                %for low frequency
                for freq = 1:length(list_low_freq{1}.freq)
                    all_session_low_freq{nvr,t,freq} =[all_session_low_freq{nvr,t,freq} list_low_freq{nvr,f}.powspctrm(chan,freq,t)];
                    all_session_low_freq_coh{nvr,t,freq} =[all_session_low_freq_coh{nvr,t,freq} list_low_freq_coh{nvr,f}.cohspctrm(chan,freq,t)];
                end
                
                
                %for high frequency
                for freq = 1:length(list_high_freq{1}.freq)
                    all_session_high_freq{nvr,t,freq} =[all_session_high_freq{nvr,t,freq} list_high_freq{nvr,f}.powspctrm(chan,freq,t)];
                    all_session_high_freq_coh{nvr,t,freq} =[all_session_high_freq_coh{nvr,t,freq} list_high_freq_coh{nvr,f}.cohspctrm(chan,freq,t)];
                end
            end
        end
    end
end

%---Calculate averages---%
all_session_average_low_freq = cell(1,2);
all_session_average_high_freq = cell(1,2);
all_session_average_low_freq_coh = cell(1,2);
all_session_average_high_freq_coh = cell(1,2);

for nvr = 1:2
    for t = 1:length(list_high_freq{1}.time)
        
        %for low frequency
        for freq = 1:length(list_low_freq{1}.freq)
            all_session_average_low_freq{nvr}(freq,t) = nanmean(all_session_low_freq{nvr,t,freq});
            all_session_average_low_freq_coh{nvr}(freq,t) = nanmean(all_session_low_freq_coh{nvr,t,freq});
        end
        
        %for high frequency
        for freq = 1:length(list_high_freq{1}.freq)
            all_session_average_high_freq{nvr}(freq,t) = nanmean(all_session_high_freq{nvr,t,freq});
            all_session_average_high_freq_coh{nvr}(freq,t) = nanmean(all_session_high_freq_coh{nvr,t,freq});
        end
    end
end


%% figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Summary Figures\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot average LFP waveforms---%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smval = [];
time = (-2048:2048)/1000;
figure
hold on
dofill(time,all_session_average_waveform{1},'blue',smval);
dofill(time,all_session_average_waveform{2},'red',smval);
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.5 0.8])
xlabel('Time From Saccade Start')
ylabel('LFP (uV)')
legend('Novel','Repeat')
title(sprintf(['List Saccade Locked LFPs for novel vs. repeat n_{electrode} = '...
    num2str(size(all_session_average_waveform{1},1))]));
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
grid(gca,'minor')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Frequency Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


freq_time = list_low_freq{1}.time;
low_freq = list_low_freq{1}.freq;
high_freq = list_high_freq{1}.freq;


yl = NaN(2,2);

figure
%---high Frequency for novel saccades---%
subplot(2,3,1)
imagesc(freq_time,high_freq,all_session_average_high_freq{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Novel')

%---high Frequency for repeat saccades---%
subplot(2,3,2)
imagesc(freq_time,high_freq,all_session_average_high_freq{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Repeat')

%---high Frequency rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:2
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(high_freq) max(high_freq)],'w--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

%---high Freq Novel-Repeat---%
subplot(2,3,3)
imagesc(freq_time,high_freq,all_session_average_high_freq{1}-all_session_average_high_freq{2})
hold on
plot([0 0],[min(high_freq) max(high_freq)],'w--')
hold off
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
title('Novel-Repeat')
colorbar, axis xy
colormap jet


yl = NaN(2,2);

%---Low Frequency for novel saccades---%
subplot(2,3,4)
imagesc(freq_time,low_freq,all_session_average_low_freq{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Novel')

%---Low Frequency for repeat saccades---%
subplot(2,3,5)
imagesc(freq_time,low_freq,all_session_average_low_freq{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Repeat')

%---Low Frequency rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 4:5
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(low_freq) max(low_freq)],'w--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

%---Low Freq Novel-Repeat---%
subplot(2,3,6)
imagesc(freq_time,low_freq,all_session_average_low_freq{1}-all_session_average_low_freq{2})
hold on
plot([0 0],[min(low_freq) max(low_freq)],'w--')
hold off
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
title('Novel-Repeat')
colorbar, axis xy
colormap jet

subtitle('List Saccades Novel vs Repeat: Power Spectral Analysis')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Coherence Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yl = NaN(2,2);

figure

%---high Frequency for novel saccades---%
subplot(2,3,1)
imagesc(freq_time,high_freq,all_session_average_high_freq_coh{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Novel')

%---high Frequency for repeat saccades---%
subplot(2,3,2)
imagesc(freq_time,high_freq,all_session_average_high_freq_coh{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Repeat')

%---high Frequency rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:2
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(high_freq) max(high_freq)],'w--')
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

%---high Freq Novel-Repeat---%
subplot(2,3,3)
imagesc(freq_time,high_freq,all_session_average_high_freq_coh{1}-all_session_average_high_freq_coh{2})
hold on
plot([0 0],[min(high_freq) max(high_freq)],'w--')
hold off
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
title('Novel-Repeat')
colorbar, axis xy
colormap jet


%---Low Frequency for novel saccades---%
subplot(2,3,4)
imagesc(freq_time,low_freq,all_session_average_low_freq_coh{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Novel')

%---Low Frequency for repeat saccades---%
subplot(2,3,5)
imagesc(freq_time,low_freq,all_session_average_low_freq_coh{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Repeat')

%---Low Frequency rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 4:5
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(low_freq) max(low_freq)],'w--')
    line([0 0],get(gca,'ylim'),'Color','w','LineWidth',2)
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

%---Low Freq Novel-Repeat---%
subplot(2,3,6)
imagesc(freq_time,low_freq,all_session_average_low_freq_coh{1}-all_session_average_low_freq_coh{2})
hold on
plot([0 0],[min(low_freq) max(low_freq)],'w--')
hold off
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
title('Novel-Repeat')
colorbar, axis xy
colormap jet


subtitle('List Saccades Novel vs Repeat: Coherence Analysis')



