% % written Seth Konig 4/29/16
% % grab analyzed data across all sessions and monkeys and get the averages
% % across all sessions and do some statistical analysis
% 
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Data\';

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

files = files(1);

sequence_average_LFP = cell(2,length(files));
sequence_low_freq = cell(2,length(files));
sequence_high_freq = cell(2,length(files));
sequence_low_freq_coh = cell(2,length(files));
sequence_high_freq_coh = cell(2,length(files));

sequence_average_LFP = cell(1,length(files));
for f = 1:length(files)
    disp([ 'Processing file ' a{files(f)}(1:8)])
    load([data_dir a{files(f)}],'sequence_saccade_aligneddata','stretched_time_to_fixation','listsq_ITI_aligneddata');

    %in case too_late was implemented too late
    stretched_time_to_fixation = stretched_time_to_fixation(1:length(sequence_saccade_aligneddata.trial));

    if strcmpi(a{files(f)}(1:2),'PW')
        predict_rt = 156;
    else %for Tobii
        predict_rt = 138;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%-Remove Trials with NaNs--%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %in case it was not implemented before (i.e. during first import)

     %---Saccades Sequence trials---%
    %remove trials with NaNs
    remove_trials = [];
    for t = 1:length(sequence_saccade_aligneddata.trial)
        if sum(sum(isnan(sequence_saccade_aligneddata.trial{t}))) >1
            remove_trials =[remove_trials t];
        end
    end
    sequence_saccade_aligneddata.time(remove_trials) = [];
    sequence_saccade_aligneddata.trial(remove_trials) = [];
    sequence_saccade_aligneddata.sampleinfo(remove_trials,:) = [];

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

    %---For Reactive Sequence Trials---%
    fixind = find(stretched_time_to_fixation < predict_rt);
    temp = sequence_saccade_aligneddata;
    temp.trial(fixind) = [];
    temp.time(fixind) = [];

    %Get Average Waveform
    timelock = ft_timelockanalysis(cfg,temp);
    sequence_average_LFP{1,f} = timelock.avg;

    %Get Low and High Frequency Power Spectra
    freq = lfp_powerspectrum(temp,'low','all',time_window);
    %normalize to baseline power by frequency and channel
    sequence_low_freq{1,f} = Apply_Seths_BaselinePower(freq,list_low_freq_ITI_BaseLinePower,'relchange');

    freq = lfp_powerspectrum(temp,'high','all',time_window);
    %normalize to baseline power by frequency and channel
    sequence_high_freq{1,f} = Apply_Seths_BaselinePower(freq,list_high_freq_ITI_BaseLinePower,'relchange');

    %Get Low and High Frequency Coherence
    sequence_low_freq_coh{1,f} = lfp_phasecoherence(temp,'low','all',time_window);
    sequence_high_freq_coh{1,f} = lfp_phasecoherence(temp,'high','all',time_window);
    
    %---For Predictive Sequence Trials---%
    fixind = find(stretched_time_to_fixation > predict_rt);
    if length(fixind) < 10 %way too few trials
       continue
    end
    temp = sequence_saccade_aligneddata;
    temp.trial(fixind) = [];
    temp.time(fixind) = [];

    %Get Average Waveform
    timelock = ft_timelockanalysis(cfg,temp);
    sequence_average_LFP{2,f} = timelock.avg;

    %Get Low and High Frequency Power Spectra
    freq = lfp_powerspectrum(temp,'low','all',time_window);
    %normalize to baseline power by frequency and channel
    sequence_low_freq{2,f} = Apply_Seths_BaselinePower(freq,list_low_freq_ITI_BaseLinePower,'relchange');

    freq = lfp_powerspectrum(temp,'high','all',time_window);
    %normalize to baseline power by frequency and channel
    sequence_high_freq{2,f} = Apply_Seths_BaselinePower(freq,list_high_freq_ITI_BaseLinePower,'relchange');

    %Get Low and High Frequency Coherence
    sequence_low_freq_coh{2,f} = lfp_phasecoherence(temp,'low','all',time_window);
    sequence_high_freq_coh{2,f} = lfp_phasecoherence(temp,'high','all',time_window);
end
save([data_dir '2ListSQ_Sequence_Reactive_vs_Predictive_saccade_aligned_averaged_across_sessions.mat'],...
    'sequence_average_LFP','sequence_low_freq',...
    'sequence_high_freq','sequence_low_freq_coh','sequence_high_freq_coh');
%
%%
%---Sequence Variables---%
all_sequence_averages_waveforms = cell(1,2);
sequence_all_session_low_freq = cell(2,length(sequence_low_freq{1}.time),length(sequence_low_freq{1}.freq));
sequence_all_session_high_freq = cell(2,length(sequence_high_freq{1}.time),length(sequence_high_freq{1}.freq));
sequence_all_session_low_freq_coh = cell(2,length(sequence_low_freq{1}.time),length(sequence_low_freq{1}.freq));
sequence_all_session_high_freq_coh = cell(2,length(sequence_high_freq{1}.time),length(sequence_high_freq{1}.freq));


for f = 1:length(sequence_average_LFP)-1 %last file doesn't necessary work properly
    for rp = 1:2 %reactive or predictive 
        all_sequence_averages_waveforms{rp} = [all_sequence_averages_waveforms{rp}; sequence_average_LFP{rp,f}];%all "waveforms"
        for chan = 1:size(sequence_low_freq{rp,f}.powspctrm,1)
            for t = 1:length(sequence_high_freq{1}.time)

                %for low frequency
                for freq = 1:length(sequence_low_freq{1}.freq)
                    sequence_all_session_low_freq{rp,t,freq} =[sequence_all_session_low_freq{rp,t,freq} sequence_low_freq{rp,f}.powspctrm(chan,freq,t)];
                    sequence_all_session_low_freq_coh{rp,t,freq} =[sequence_all_session_low_freq_coh{rp,t,freq} sequence_low_freq_coh{rp,f}.cohspctrm(chan,freq,t)];
                end

                %for high frequency
                for freq = 1:length(sequence_high_freq{1}.freq)
                    sequence_all_session_high_freq{rp,t,freq} =[sequence_all_session_high_freq{rp,t,freq} sequence_high_freq{rp,f}.powspctrm(chan,freq,t)];
                    sequence_all_session_high_freq_coh{rp,t,freq} =[sequence_all_session_high_freq_coh{rp,t,freq} sequence_high_freq_coh{rp,f}.cohspctrm(chan,freq,t)];
                end
            end
        end
    end
end

%---Calculate averages---%
sequence_all_session_average_low_freq = cell(1,2);
sequence_all_session_average_high_freq = cell(1,2);
sequence_all_session_average_low_freq_coh = cell(1,2);
sequence_all_session_average_high_freq_coh = cell(1,2);

%for list
for rp = 1:2
    for t = 1:length(sequence_high_freq{1}.time)

        %for low frequency
        for freq = 1:length(sequence_low_freq{1}.freq)
            sequence_all_session_average_low_freq{rp}(freq,t) = nanmean(sequence_all_session_low_freq{rp,t,freq});
            sequence_all_session_average_low_freq_coh{rp}(freq,t) = nanmean(sequence_all_session_low_freq_coh{rp,t,freq});
        end

        %for high frequency
        for freq = 1:length(sequence_high_freq{1}.freq)
            sequence_all_session_average_high_freq{rp}(freq,t) = nanmean(sequence_all_session_high_freq{rp,t,freq});
            sequence_all_session_average_high_freq_coh{rp}(freq,t) = nanmean(sequence_all_session_high_freq_coh{rp,t,freq});
        end
    end
end

%%
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\Summary Figures\';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot average LFP waveforms---%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = (-2048:2048)/1000;
smval = [];

figure
hold on
dofill(time,all_sequence_averages_waveforms{1},'red',smval);
dofill(time,all_sequence_averages_waveforms{2},'blue',smval);
yl = ylim;
plot([0 0],[yl(1) yl(2)],'k--')
hold off
xlim([-0.5 0.8])
xlabel('Time From Saccade Start')
ylabel('LFP (uV)')
legend('Reactive','Predictive')
title(sprintf(['Saccade Locked LFPs For Sequences n_{electrode} = '...
    num2str(size(all_sequence_averages_waveforms{1},1))]));
grid on
set(gca,'XMinorTick','on','YMinorTick','on')
grid on
grid(gca,'minor')

save_and_close_fig(figure_dir,'Sequence_Predictive_vs_Reactive-Waveforms')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Frequency Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freq_time = sequence_low_freq{1}.time;
low_freq = sequence_low_freq{1}.freq;
high_freq = sequence_high_freq{1}.freq;


yl = NaN(2,2);

figure
%---Reactive high Frequency Saccades---%
subplot(2,3,1)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Reactive')

%---Predictive high Frequency Saccades---%
subplot(2,3,2)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Predictive')

%---Reactive-Predictive---%
subplot(2,3,3)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq{1}-sequence_all_session_average_high_freq{2})
hold on
plot([0 0],[min(high_freq) max(high_freq)],'k--')
hold off
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Reactive-Predictive')
colorbar, axis xy
colormap jet

%---rescale and add line plot---%
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

yl = NaN(2,2);
%---Reactive Low Frequency Saccades---%
subplot(2,3,4)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Reactive')

%---Predictive Low Frequency Saccades---%
subplot(2,3,5)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Predictive')

%---Reactive-Predictives---%
subplot(2,3,6)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq{1}-sequence_all_session_average_low_freq{2})
hold on
plot([0 0],[min(low_freq) max(low_freq)],'k--')
hold off
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Reactive-Predictive')
colorbar, axis xy
colormap jet

%---rescale and add line plot---%
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

save_and_close_fig(figure_dir,'Sequence_Predictive_vs_Reactive-Frequency Power Spectrum')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---Plot Coherence Analysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yl = NaN(2,2);

figure

%---High Freqency Reactive Saccades---%
subplot(2,3,1)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq_coh{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Reactive')

%---High Freqency Predictive Saccades---%
subplot(2,3,2)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq_coh{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Predictive')

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 1:2
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(high_freq) max(high_freq)],'k--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

%---High Freqency Reactive-Predictive Saccades---%
subplot(2,3,3)
imagesc(freq_time,high_freq,sequence_all_session_average_high_freq_coh{1}-sequence_all_session_average_high_freq_coh{2})
hold on
plot([0 0],[min(high_freq) max(high_freq)],'k--')
hold off
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
title('Reactive-Predictive')
colorbar, axis xy
colormap jet


yl = NaN(2,2);
%---low Freqency Reactive Saccades---%
subplot(2,3,4)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq_coh{1})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(1,:) = caxis;
title('Reactive')

%---low Freqency Predictive Saccades---%
subplot(2,3,5)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq_coh{2})
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Predictive')

%---rescale and add line plot---%
ymin = min(yl(:,1));
ymax = max(yl(:,2));
for sb = 4:5
    subplot(2,3,sb)
    hold on
    plot([0 0],[min(low_freq) max(low_freq)],'k--')
    colorbar, axis xy
    colormap jet
    hold off
    caxis([ymin ymax])
end

%---low Freqency Reactive-Predictive Saccades---%
subplot(2,3,6)
imagesc(freq_time,low_freq,sequence_all_session_average_low_freq_coh{1}-sequence_all_session_average_low_freq_coh{2})
hold on
plot([0 0],[min(low_freq) max(low_freq)],'k--')
hold off
ylabel('Frequency (Hz)')
xlim([-0.75 0.75])
xlabel('Time from Saccade (sec)')
yl(2,:) = caxis;
title('Reactive-Predictive')
colorbar, axis xy
colormap jet

save_and_close_fig(figure_dir,'Sequence_Predictive_vs_Reactive-Coherence Spectrum')

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

predictive_freq_pre_post_ps = NaN(1,size(band,1));
predictive_freq_means = NaN(2,size(band,1));
predictive_freq_stds = NaN(2,size(band,1));

reactive_freq_pre_post_ps = NaN(1,size(band,1));
reactive_freq_means = NaN(2,size(band,1));
reactive_freq_stds = NaN(2,size(band,1));


colors = {'red','blue'};

%for frequency spectrum analysis
figure
for b = 1:size(band,1);
    
    %---For List---%
    if b <= 3 %for the low frequencies
        freq_ind = find(low_freq >= band(b,1) & low_freq <= band(b,2));
        
        subplot(2,3,b)
        hold on
        for rp = 1:2
            %grab the data
            temp = cell(t,length(freq_ind));
            for f = 1:length(freq_ind)
                for t = 1:size(sequence_all_session_low_freq,2)
                    temp{t,f} = sequence_all_session_low_freq{rp,t,freq_ind(f)};
                end
            end
            %average across the frequencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(sequence_all_session_low_freq,2));
            for t = 1:size(sequence_all_session_low_freq,2)
                temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
            end
            dofill(freq_time,temp2,colors{rp},smval);
            
            pre_window = nanmean(temp2(:,pre_time),2);
            post_window = nanmean(temp2(:,post_time));
            
            if rp == 1
                reactive_freq_means(2,b) = nanmean(post_window);
                reactive_freq_stds(2,b) = nanstd(post_window);
                reactive_freq_means(1,b) = nanmean(pre_window);
                reactive_freq_stds(1,b) = nanstd(pre_window);
                [~,reactive_freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
            else
                predictive_freq_means(2,b) = nanmean(post_window);
                predictive_freq_stds(2,b) = nanstd(post_window);
                predictive_freq_means(1,b) = nanmean(pre_window);
                predictive_freq_stds(1,b) = nanstd(pre_window);
                [~,predictive_freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
            end
        end

    else %for the high frequencies
        freq_ind = find(high_freq >= band(b,1) & high_freq <= band(b,2));

        subplot(2,3,b)
        hold on
        for rp = 1:2
            %grab the data
            temp = cell(t,length(freq_ind));
            for f = 1:length(freq_ind)
                for t = 1:size(sequence_all_session_high_freq,2)
                    temp{t,f} = sequence_all_session_high_freq{rp,t,freq_ind(f)};
                end
            end
            %average across the frequencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(sequence_all_session_high_freq,2));
            for t = 1:size(sequence_all_session_high_freq,2)
               temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_ind))')';
            end
            dofill(freq_time,temp2,colors{rp},smval); 
            
              pre_window = nanmean(temp2(:,pre_time),2);
            post_window = nanmean(temp2(:,post_time));
            
            if rp == 1
                reactive_freq_means(2,b) = nanmean(post_window);
                reactive_freq_stds(2,b) = nanstd(post_window);
                reactive_freq_means(1,b) = nanmean(pre_window);
                reactive_freq_stds(1,b) = nanstd(pre_window);
                [~,reactive_freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
            else
                predictive_freq_means(2,b) = nanmean(post_window);
                predictive_freq_stds(2,b) = nanstd(post_window);
                predictive_freq_means(1,b) = nanmean(pre_window);
                predictive_freq_stds(1,b) = nanstd(pre_window);
                [~,predictive_freq_pre_post_ps(b)] = ttest2(post_window,pre_window);
            end
            

        end
    end

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
       legend('Reactive','Predictive') 
    end
end

%bar graph for power before and after saccade for reactive
subplot(4,3,9) %for list only 500-600 ms
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[reactive_freq_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[reactive_freq_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[reactive_freq_means(1,b) reactive_freq_means(2,b)],[reactive_freq_stds(1,b) reactive_freq_stds(2,b)]./nsamp)
    if reactive_freq_pre_post_ps(b) < 0.05
       plot(b,0.2,'*k') 
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Relative Power')
title('Reactive')

%bar graph for power before and after saccade for predictive
subplot(4,3,12)
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[predictive_freq_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[predictive_freq_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[predictive_freq_means(1,b) predictive_freq_means(2,b)],[predictive_freq_stds(1,b) predictive_freq_stds(2,b)]./nsamp)
    if predictive_freq_pre_post_ps(b) < 0.25
       plot(b,0.5,'*k') 
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Relative Power')
title('Predictive')

save_and_close_fig(figure_dir,'Sequence_Predictive_vs_Reactive-Line Plot Frequency Power Spectrum')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%---freq_cohuency Statistical Anlaysis---%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
twin = .150;
pre_time = find(freq_time >= -twin & freq_time  < 0);
post_time = find(freq_time >= 0 & freq_time < twin);
band  = [4 12; %theta
    12 16; %alpha
    16 30; %beta
    30 50; %low gamma
    50 120]; %high gamma

predictive_freq_coh_pre_post_ps = NaN(1,size(band,1));
predictive_freq_coh_means = NaN(2,size(band,1));
predictive_freq_coh_stds = NaN(2,size(band,1));

reactive_freq_coh_pre_post_ps = NaN(1,size(band,1));
reactive_freq_coh_means = NaN(2,size(band,1));
reactive_freq_coh_stds = NaN(2,size(band,1));


colors = {'red','blue'};

%for freq_cohuency spectrum analysis
figure
for b = 1:size(band,1);
    
    %---For List---%
    if b <= 3 %for the low freq_cohuencies
        freq_coh_ind = find(low_freq >= band(b,1) & low_freq <= band(b,2));
        
        subplot(2,3,b)
        hold on
        for rp = 1:2
            %grab the data
            temp = cell(t,length(freq_coh_ind));
            for f = 1:length(freq_coh_ind)
                for t = 1:size(sequence_all_session_low_freq_coh,2)
                    temp{t,f} = sequence_all_session_low_freq_coh{rp,t,freq_coh_ind(f)};
                end
            end
            %average across the freq_cohuencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(sequence_all_session_low_freq_coh,2));
            for t = 1:size(sequence_all_session_low_freq_coh,2)
                temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_coh_ind))')';
            end
            dofill(freq_time,temp2-min(nanmean(temp2)),colors{rp},smval);
            
            pre_window = nanmean(temp2(:,pre_time),2);
            post_window = nanmean(temp2(:,post_time));
            
            if rp == 1
                reactive_freq_coh_means(2,b) = nanmean(post_window);
                reactive_freq_coh_stds(2,b) = nanstd(post_window);
                reactive_freq_coh_means(1,b) = nanmean(pre_window);
                reactive_freq_coh_stds(1,b) = nanstd(pre_window);
                [~,reactive_freq_coh_pre_post_ps(b)] = ttest2(post_window,pre_window);
            else
                predictive_freq_coh_means(2,b) = nanmean(post_window);
                predictive_freq_coh_stds(2,b) = nanstd(post_window);
                predictive_freq_coh_means(1,b) = nanmean(pre_window);
                predictive_freq_coh_stds(1,b) = nanstd(pre_window);
                [~,predictive_freq_coh_pre_post_ps(b)] = ttest2(post_window,pre_window);
            end
        end

    else %for the high freq_cohuencies
        freq_coh_ind = find(high_freq >= band(b,1) & high_freq <= band(b,2));

        subplot(2,3,b)
        hold on
        for rp = 1:2
            %grab the data
            temp = cell(t,length(freq_coh_ind));
            for f = 1:length(freq_coh_ind)
                for t = 1:size(sequence_all_session_high_freq_coh,2)
                    temp{t,f} = sequence_all_session_high_freq_coh{rp,t,freq_coh_ind(f)};
                end
            end
            %average across the freq_cohuencies in the band but keep
            %electrode-session seperate
            temp2 = NaN(size(temp{1},2),size(sequence_all_session_high_freq_coh,2));
            for t = 1:size(sequence_all_session_high_freq_coh,2)
               temp2(:,t) = nanmean(reshape(cell2mat(temp(t,:)),size(temp{1},2),length(freq_coh_ind))')';
            end
            dofill(freq_time,temp2-min(nanmean(temp2)),colors{rp},smval); 
            
            pre_window = nanmean(temp2(:,pre_time),2);
            post_window = nanmean(temp2(:,post_time));
            
            if rp == 1
                reactive_freq_coh_means(2,b) = nanmean(post_window);
                reactive_freq_coh_stds(2,b) = nanstd(post_window);
                reactive_freq_coh_means(1,b) = nanmean(pre_window);
                reactive_freq_coh_stds(1,b) = nanstd(pre_window);
                [~,reactive_freq_coh_pre_post_ps(b)] = ttest2(post_window,pre_window);
            else
                predictive_freq_coh_means(2,b) = nanmean(post_window);
                predictive_freq_coh_stds(2,b) = nanstd(post_window);
                predictive_freq_coh_means(1,b) = nanmean(pre_window);
                predictive_freq_coh_stds(1,b) = nanstd(pre_window);
                [~,predictive_freq_coh_pre_post_ps(b)] = ttest2(post_window,pre_window);
            end
            

        end
    end

end

band_str = {'Theta','Alpha','Beta','Low Gamma','High Gamma'};
nsamp = sqrt(size(temp2,1));
for sb = 1:5
    subplot(2,3,sb)
    yl = ylim;
    if yl(1) < 0
        yl(1) = 0;
        ylim(yl)
    end
    plot([0 0],[yl(1) yl(2)],'k--')
    hold off
    xlim([-0.75 0.75])
    xlabel('Time From Saccade Start')
    ylabel('Baseline-Corrected Coherence')
    title([band_str{sb} ': ' num2str(band(sb,1)) '-'  num2str(band(sb,2)) ' Hz'])
    
    if sb == 5
       legend('Reactive','Predictive') 
    end
end

%bar graph for power before and after saccade for reactive
subplot(4,3,9) %for list only 500-600 ms
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[reactive_freq_coh_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[reactive_freq_coh_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[reactive_freq_coh_means(1,b) reactive_freq_coh_means(2,b)],[reactive_freq_coh_stds(1,b) reactive_freq_coh_stds(2,b)]./nsamp)
    if reactive_freq_coh_pre_post_ps(b) < 0.05
       plot(b,0.2,'*k') 
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Baseline-Corrected Coherence')
title('Reactive')

%bar graph for power before and after saccade for predictive
subplot(4,3,12)
hold on
for b = 1:length(band)
    bp = bar([b-0.25],[predictive_freq_coh_means(1,b)],0.3); %pre
    set(bp,'FaceColor','r','EdgeColor','k')
    bp = bar([b+0.25],[predictive_freq_coh_means(2,b)],0.3); %post
    set(bp,'FaceColor','g','EdgeColor','k')
    errorb([b-0.25 b+0.25],[predictive_freq_coh_means(1,b) predictive_freq_coh_means(2,b)],[predictive_freq_coh_stds(1,b) predictive_freq_coh_stds(2,b)]./nsamp)
    if predictive_freq_coh_pre_post_ps(b) < 0.25
       plot(b,0.5,'*k') 
    end
end
set(gca,'Xtick',1:b)
set(gca,'XtickLabel',band_str)
ylabel('Baseline-Corrected Coherence')
title('Predictive')

%save_and_close_fig(figure_dir,'Absolute_Sequence_Predictive_vs_Reactive-Line Plot Coherence Spectrum')
save_and_close_fig(figure_dir,'Corrected_Sequence_Predictive_vs_Reactive-Line Plot Coherence Spectrum')

