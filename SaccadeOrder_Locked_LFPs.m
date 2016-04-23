% Written to determine if LFPs are obviously locked to saccades
% Written November 21, 2014 by Seth  Konig
% modified code 10/26/15 to look at differneces in saccade order


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Analysis for Images%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
% fixation_durations = {};
%
% cd(dir);
% a = what;
% m = a.mat;
% set = 1;
% average_LFP = [];
% for file = 1:size(m,1);
%     average_LFP{set} = cell(35,4); %set by electrode by sequence vs image
%     if ~isempty(strfind(m{file},'preprocessed'))
%         fixationstats = [];
%         load(m{file},'fixationstats','cfg','item_set','data','hdr');
%         if ~isempty(strfind(item_set,'ckwenzr'))
%             continue;
%         end
%         lfpchans = find_desired_channels(hdr,data,'LFP');
%         if ~isempty(fixationstats) %so is ListsQ
%             [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
%             [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
%             num_sacs = 0;
%             for t = 1:length(fixationstats);
%                 if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
%                     trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
%                     imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start;%onlg get saccades while image is up
%                     imgoff =  cfg.trl(t).alltim(cfg.trl(t).allval == 24)-trial_start;%
%                     dataend =  cfg.trl(t).alltim(cfg.trl(t).allval == 20)-trial_start;%to make sure I have enough data for 500 ms after saccade
%                     saccadetimes = fixationstats{t}.saccadetimes;
%                     if ~isempty(saccadetimes)
%                         tooearly = find(saccadetimes(1,:) <= imgon+500);%also want to ignore the 1st saccade since evoked by image
%                         saccadetimes(:,tooearly) = [];
%                         toolate = find(saccadetimes(1,:) >= imgoff);
%                         saccadetimes(:,toolate) = [];
%                         tootoolate =  find(saccadetimes(1,:) >= dataend-500);
%                         saccadetimes(:,tootoolate) = [];
%                         maxsacs = size(saccadetimes,2);
%                         maxsacs(maxsacs > 35) = 35;
%                         for s = 1:maxsacs
%                             for l = 1:4
%                                 average_LFP{set}{s,l} = [average_LFP{set}{s,l}; data(lfpchans(l)).values{t}(saccadetimes(1,s):saccadetimes(1,s)+499)];
%                             end
%                         end
%                         num_sacs = num_sacs + maxsacs;
%                     end
%                 end
%             end
%             set = set+1;
%         end
%     end
% end
% %%
% % Normalize LFPs on each channel by RMS
% setRMS = cell(1,length(average_LFP));
% for set = 1:length(average_LFP)
%     RMS = NaN(4,size(average_LFP{set},1));
%     for l = 1:4
%         for s = 1:size(average_LFP{set},1)
%             RMS(l,s) = sqrt(sum(sum(average_LFP{set}{s,l}.^2))/numel(average_LFP{set}{s,l}));
%         end
%     end
%     setRMS{set} = nanmean(RMS');
% end
%
% for set = 1:length(average_LFP)
%     figure
%     for l = 1:4
%         subplot(2,2,l)
%         hold all
%         for s = 1:size(average_LFP{set},1)
%             plot(mean(average_LFP{set}{s,l}/setRMS{set}(l)))
%         end
%     end
% end
%
% sacnum = [];
% enum = [];
% all_LFP = [];
% for set = 1:length(average_LFP)
%     for l = 1:4
%         for s = 1:size(average_LFP{set},1)
%             all_LFP = [all_LFP; average_LFP{set}{s,l}/setRMS{set}(l)];
%             sacnum = [sacnum,s*ones(1,size(average_LFP{set}{s,l},1))];
%             enum =  [enum,l*ones(1,size(average_LFP{set}{s,l},1))];
%         end
%     end
% end
% [pc,score,latent,tsquare] = princomp(all_LFP);
%
% [p,~,STATS] = anovan(score(:,1),{sacnum',enum'});
% multcompare(STATS)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Analysis for Sequence%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';

twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
smval = 120;%gaussian 1/2 width for smoothing
predicted_thresh = 10;% percent of saccades that must be < 150 ms to item to constitue as a learned-predicted item
% for Vivian 5% was the false positive rate for 150 ms so wanted at least double the False positive rate
predicted_rt = 150;%maximum "reaction time" for what constitutes as predictive, everything else is reactive
numshuffs = 100; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;

cd(data_dir);
a = what;
m = a.mat;
set = 1;
average_LFP = [];
for file = 1:size(m,1);
    average_LFP{set} = cell(4,4); %set by electrode by sequence vs image
    if ~isempty(strfind(m{file},'preprocessed'))
        fixationstats = [];
        load(m{file},'fixationstats','cfg','item_set','data','hdr');
        lfpchans = length(data)-5:length(data)-2;
        if ~isempty(fixationstats) %so is ListsQ
            
            %get important task specific information
            [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
            
            %preallocate space and parallel structure of cfg
            successful_sequence_trials = NaN(1,length(cfg.trl));
            which_sequence = NaN(1,length(cfg.trl));
            for t = 1:length(cfg.trl);
                if sum(cfg.trl(t).allval == 3) >= 6; %in which sequence trials were rewarded
                    which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
                    successful_sequence_trials(t) = t;
                end
            end
            successful_sequence_trials = laundry(successful_sequence_trials);
            which_sequence = laundry(which_sequence);
            
            num_trials = length(successful_sequence_trials);
            event_times = NaN(length(successful_sequence_trials),8);
            startindex = NaN(1,length(successful_sequence_trials)); %index for all data of trial start
            for t = 1:num_trials
                trial_start = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == trial_start_code);
                startindex(t) = trial_start;
                for event = 1:length(event_codes);
                    event_times(t,event) = cfg.trl(successful_sequence_trials(t)).alltim(cfg.trl(successful_sequence_trials(t)).allval == event_codes(event))-trial_start;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---process eye data locked to trial events---%%%
            fixationstats = fixationstats(successful_sequence_trials);
            cfg.trl = cfg.trl(successful_sequence_trials);
            
            saccade_start_time = NaN(length(fixationstats),4);%when did saccade to item start
            fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
            time_to_fixation = NaN(length(fixationstats),4); %how fast did they get to it the first time around
            for trial = 1:num_trials
                locs = sequence_locations{which_sequence(trial)};
                
                %convert to DVA for this analysis
                locs(1,:) = (locs(1,:)-400)/24;
                locs(2,:) = (locs(2,:)-300)/24;
                fixationstats{trial}.fixations(1,:) = (fixationstats{trial}.fixations(1,:)-400)/24;
                fixationstats{trial}.fixations(2,:) = (fixationstats{trial}.fixations(2,:)-300)/24;
                
                event_codes = cfg.trl(trial).allval;
                event_codes(event_codes == 100)= 0;
                event_codes(1) = 100;%eye data starts for recording right away
                event_times = cfg.trl(trial).alltim;
                
                trialdata = analyze_sequence_trial(fixationstats{trial},locs,fixwin,...
                    event_codes,event_times);
                
                time_to_fixation(trial,:) = trialdata.t2f;

                fixation_numbers = trialdata.fixationnums; %fixation number for each item
                fixationtimes = fixationstats{trial}.fixationtimes;
                saccadetimes = fixationstats{trial}.saccadetimes;
                if saccadetimes(1,1) < fixationtimes(1,1); %then started with a saccade
                    sacind = 0;
                else%started with a fixation
                    sacind = -1;
                end
                for item = 1:4
                    if ~isnan(fixation_numbers(item))
                        fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
                        if fixation_numbers(item)+sacind  >= 1
                            try
                                saccade_start_time(trial,item) = saccadetimes(1,fixation_numbers(item)+sacind);
                            catch
                                saccade_start_time(trial,item) = NaN;%should only occur if saccade was off screen or a blink occured otherwise should work
                                disp('Saccade Start Time not found')
                            end
                        end
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%---Analyzing Predictive vs Reactive LFP Data---%%%
            
            %%%---First Remove line noise & it's harmonics---%%%
            Fline       = 60; %removes harmonics so 120, 180, etc.
            LFPchannels = find_desired_channels(hdr,data,'LFP');
            Fs = hdr.Fs;
            
            sequence_LFPdata = cell(4,num_trials);
            for trial = 1:num_trials
                
                LFPdata = [];
                for channel = 1:length(LFPchannels)
                    LFPdata = [LFPdata; data(LFPchannels(channel)).values{successful_sequence_trials(trial)}];
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
                    sequence_LFPdata{channel,trial} = LFPdata(channel,:);
                end
            end
            
            %%%---Align LFPs to Saccades and Fixations---%%%
            %channel data alternate every 4th row
            fixation_aligned_LFP = cell(1,4);
            saccade_aligned_LFP = cell(1,4);
            predicted_vs_reactive = NaN(4,4*num_trials);
            sequence = NaN(4,4*num_trials);
            for c = 1:4
                for trial = 1:num_trials
                    %fixt = fixation_start_time(trial,c);
                    sact = saccade_start_time(trial,c);
                    for channel = 1:4
                        %                         if fixt > 512
                        %                             fixation_aligned_LFP{c}(4*(trial-1)+channel,:) = sequence_LFPdata{channel,trial}(fixt-512:fixt+511);
                        %                         end
                        if sact > 512 && sact+511 <= length(sequence_LFPdata{channel,trial})
                            average_LFP{set}{c,channel} = [average_LFP{set}{c,channel}; sequence_LFPdata{channel,trial}(sact-512:sact+511)];
                        end
                        %predicted_vs_reactive(c,4*(trial-1)+channel) = predicted(trial,c);
                        %sequence(c,4*(trial-1)+channel) = which_sequence(trial);
                    end
                end
            end
            set = set +1;
        end
    end
end
%%

% Normalize LFPs on each channel by RMS
setRMS = cell(1,length(average_LFP));
for set = 1:length(average_LFP)
    RMS = NaN(4,size(average_LFP{set},1));
    for l = 1:4
        for s = 1:size(average_LFP{set},1)
            RMS(l,s) = sqrt(sum(sum(average_LFP{set}{s,l}.^2))/numel(average_LFP{set}{s,l}));
        end
    end
    setRMS{set} = nanmean(RMS');
end
%%
for sets = 1:length(average_LFP)
    figure
    for l = 1:4
        subplot(2,2,l)
        hold all
        for s = 1:size(average_LFP{sets},1)
            plot(mean(average_LFP{sets}{s,l}/setRMS{sets}(l)))
        end
        plot([513 513],[-0.75 0.75],'k--')
        legend('1','2','3','4')
        xlabel('Saccade Onset (ms)')
        set(gca,'Xtick',[0:128:1024])
        set(gca,'XtickLabel',num2cell([0:128:1024]-512))
        xlim([0 1024])
    end
end
%%
sacnum = [];
enum = [];
all_LFP = [];
for set = 1:length(average_LFP)
    for l = 1:4
        for s = 1:size(average_LFP{set},1)
            all_LFP = [all_LFP; average_LFP{set}{s,l}/setRMS{set}(l)];
            sacnum = [sacnum,s*ones(1,size(average_LFP{set}{s,l},1))];
            enum =  [enum,l*ones(1,size(average_LFP{set}{s,l},1))];
        end
    end
end
[pc,score,latent,tsquare] = princomp(all_LFP);
%%
[p,~,STATS] = anovan(score(:,2),{sacnum',enum'});
multcompare(STATS)
%%
plot3(enum,score(:,1),score(:,2),'.')
%%
pcas = score(:,1:3);
thresh = min(pcas(1:end)):(max(pcas(1:end))-min(pcas(1:end)))/100:max(pcas(1:end));
allROCs = cell(1,length(average_LFP));
figure
hold all
for set = 1
    ROCs = cell(4,4);
    for s = 1:4
        for ss = 1:size(average_LFP{set},1)
            if s > ss
                PCs = pcas(sacnum == s,1);
                PCss = pcas(sacnum == ss,1);
                
                TP = NaN(1,length(thresh)); %True positive
                FA = NaN(1,length(thresh)); %False alarm
                for ii = 1:length(thresh)
                    TP(ii) = sum(PCs > thresh(ii))/length(PCs);
                    FA(ii) = sum(PCss > thresh(ii))/length(PCss);
                end
                ROCs{s,ss} = [TP;FA];
                plot(FA,TP)
            end
        end
        allROCs{set} = ROCs;
    end
end
plot([0 1],[0 1],'k--')
hold off
axis square

%%

figure
for PC = 1:3
    p10 = prctile(pcas(:,PC),10);
    p90 = prctile(pcas(:,PC),90);
    p10i = find(pcas(:,PC) < p10);
    p90i = find(pcas(:,PC) > p90);
    
    subplot(1,3,PC)
    hold all
    plot(nanmean(all_LFP(p10i,:)))
    plot(nanmean(all_LFP(p90i,:)))
    hold off
end

%% FFT

Fs = 1000;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = 512;             % Length of signal
t = (0:L-1)*T;        % Time vector

allfft = [];

for sampl = 1:size(all_LFP,1);
    y = fft(all_LFP(sampl,:));
    y = abs(y/L);
    y = y(1:L/2+1);
    y(2:end-1) = 2*y(2:end-1);
    allfft = [allfft; y];
end

f = Fs*(0:(L/2))/L;
plot(f,mean(allfft))

%% spectrogram

allspec = zeros(65,191);
for sampl = 1:size(all_LFP,1);
    s = spectrogram(all_LFP(sampl,:),16,8,128,1e3);
    allspec = allspec+s;
end

figure
imagesc(abs(allspec(end:-1:1,:)))
ylabel('Frequency')
xlabel('Time')
