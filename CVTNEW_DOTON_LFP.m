%draft analysis looking at spikes aligned to LFPs during covert attention
%task


clar
task = 'cvtnew';

%Cortex codes
trial_start_code = 15;
cross_on_code = 35;
fixation_code = 8;
dot_on_code = 25;
dot_clrchng_code = 27;
response_code = 4;
reward_code = 3;

min_blks = 1;

twin1 = 1000;
twin2 = 2000;

%spike times vs spike phase
num_frequencies = 32;

all_LFPs = [];
all_LFP_phases = cell(1,num_frequencies);
all_mrls = cell(1,num_frequencies);

phase_time = cell(num_frequencies,300);%by unit
unit_index = 0;
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';

        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])

        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML

    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';

        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas

        %         listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end

    for sess = 1:length(session_data)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import task and unit data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [task_file,~,~,multiunit,unit_names,unit_confidence,...
            sorting_quality,~,lfp_quality,~] = get_task_data(session_data{sess},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr');

        %grap LPF channel info data
        num_trials = length(cfg.trl); %number of trials

        LFPchannels = NaN(1,4);
        try
            bad_channels = [];
            for channel = 1:4
                if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                    if  lfp_quality(channel) == 0; %if it is bad
                        bad_channels = [bad_channels channel];
                    else
                        temp = strfind(cfg.channel,['AD0' num2str(channel)]);
                        temp(cellfun(@isempty,temp)) = {0};
                        LFPchannels(channel) = find(cell2mat(temp));
                    end
                end
            end
            LFPchannels(bad_channels) = NaN;
        catch
            continue
        end

        if isempty(LFPchannels)
            continue
        end

        %grab unit data
        [multiunit,unit_stats,num_units] = get_unit_names(cfg,hdr,data,unit_names,...
            multiunit,unit_confidence,sorting_quality);
        clear unit_names
        if num_units == 0
            continue %if no units exit function
        end

        %get trials in which spikes are valid i.e. trials for which neuron is
        %stable for
        valid_trials = determine_valid_trials(task_file,valid_trials,cfg,num_units,min_blks,task);
        if all(isnan(valid_trials(1,:)))
            disp([task_file(1:8) ': no valid units could be found. Exiting function...'])
            continue %since no valid units skip analysis
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---import & reformat data so that spikes are locked to events---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([task_file(1:8) ' Aligning spike times to trial events'])

        num_trials = length(cfg.trl);
        trial_LFP = cell(1,length(LFPchannels));
        trial_phase = cell(1,length(LFPchannels));
        for chan = 1:length(LFPchannels);
            if ~isnan(LFPchannels(chan))
                trial_LFP{chan} = NaN(num_trials,twin1*2+twin2);
                for freq = 1:num_frequencies
                    trial_phase{chan,freq} = NaN(num_trials,twin1*2+twin2);
                end
            end
        end
        trial_dur = NaN(1,num_trials);

        %find which channel unit was recorded on
        LFP_unit_chan = NaN(1,num_units);
        for unit = 1:num_units
            for freq = 1:num_frequencies
                phase_time{freq,unit+unit_index} = NaN(2,10000);
            end
            LFP_unit_chan(unit) = str2double(unit_stats{1,unit}(6));
        end

        %don't analyze channels without units
        for chan = 1:4
            if sum(LFP_unit_chan == chan) == 0
                LFPchannels(chan) = NaN;
            end
        end

        index = ones(num_frequencies,num_units);
        for t = 1:num_trials
            if any(cfg.trl(t).allval == reward_code); %in which trial was successful

                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
                trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after reward
                doton = cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %when image turns on
                dotchange = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start; %when dot changes color
                monkey_response =  cfg.trl(t).alltim(cfg.trl(t).allval == response_code)-trial_start; %when release bar
                reward = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start; %reward pulses

                %                 if dotchange-doton > 2000
                %                     dotchange = doton+1999;
                %                 end
                trial_dur(t) = dotchange-doton;

                for chan = 1:length(LFPchannels)
                    if ~isnan(LFPchannels(chan))
                        LFPs = data(LFPchannels(chan)).values{t}; %spike trains for this trial
                        LFPs = LFPs(doton-twin1:dotchange+twin1);
                        trial_LFP{chan}(t,1:twin1+trial_dur(t)) =  LFPs(1:twin1+trial_dur(t));


                        [~,trialphase,wfq] = waveletanalysis(LFPs,'n_fq',num_frequencies);
                        for freq = 1:num_frequencies
                            trial_phase{chan,freq}(t,1:twin1+trial_dur(t)) = trialphase(freq,1:twin1+trial_dur(t));
                        end
                    end
                end

                for unit = 1:num_units
                    chan = LFP_unit_chan(unit);

                     if ~isnan(LFPchannels(chan))
                        for freq = 1:num_frequencies
                            spikes = find(data(unit).values{t}); %spike trains for this trial
                            spikes(spikes <= doton-twin1) = [];
                            spikes(spikes > dotchange) = [];
                            spikes =spikes-doton+twin1; %zero to time within LFP

                            for s = 1:length(spikes);
                                phase_time{freq,unit+unit_index}(:,index(freq,unit)) = [spikes(s)-doton trial_phase{chan,freq}(t,spikes(s))]; %time, phase
                                index(freq,unit) = index(freq,unit)+1;
                            end
                        end
                    end
                end


            end
        end
        unit_index = unit_index+num_units;

        for chan = 1:length(LFPchannels)
            if ~isnan(LFPchannels(chan))
                for freq = 1:num_frequencies
                    mrl = NaN(1,twin1+twin2);
                    for t=1:twin1+twin2
                        vec = trial_phase{chan,freq}(:,t);
                        vec(isnan(vec)) = [];
                        if ~isempty(vec)
                            mrl(t) = circ_r(vec);
                        end
                    end
                    all_mrls{freq} = [all_mrls{freq}; mrl];

                    all_LFP_phases{freq} = [all_LFP_phases{freq}; nanmean(trial_phase{chan,freq})];
                end

                all_LFPs = [all_LFPs; nanmean(trial_LFP{chan})];
            end
        end
    end
end
%
tm = -twin1:twin1+twin2-1;
figure
plot(tm,nanmean(all_LFPs));
yl = ylim;
hold on
plot([-300 -300],[yl(1) yl(2)],'r--')
plot([0 0],[yl(1) yl(2)],'k--')
plot([700 700],[yl(1) yl(2)],'g--')
hold off
xlabel('Time from Dot On (ms)')
ylabel('Session Average LFP')
legend('LFP','Fixation','Dot On','Min Trial Duration')
box off
xlim([-twin1 twin2])

%%
fc = 4;
w = 2*pi*fc;
fn = 1000/2;
[b a] = butter(8,w/fn,'high');

filtered_LFP = NaN(size(all_LFPs));
for l = 1:size(all_LFPs,1)
    filtered_LFP(l,:) = filtfilt(b,a,all_LFPs(l,:));
end
%%


fc = 16; % Cut off frequency
fs = 1000; % Sampling rate

[b,a] = butter(6,fc/(fs/2),'high'); % Butterworth filter of order 6
filtered_LFP = NaN(size(all_LFPs));
for l = 1:size(all_LFPs,1)
    filtered_LFP(l,1:3000) = filtfilt(b,a,all_LFPs(l,1:3000));
end

%%
%%
figure
for freq = 2:2:num_frequencies
    subplot(4,4,freq/2)
    plot(tm,nanmean(all_LFP_phases{freq}))
    xlim([-500 twin2])
    title(num2str(wfq(freq),3))
    yl = ylim;
    hold on
    plot([-300 -300],[yl(1) yl(2)],'r--')
    plot([0 0],[yl(1) yl(2)],'k--')
    plot([700 700],[yl(1) yl(2)],'g--')
    hold off
    box off
end
subtitle('Mean Phase')
%%
figure
for freq =  2:2:num_frequencies
    subplot(4,4,freq/2)
    plot(-twin1:twin2-1,nanmean(all_mrls{freq}))
    xlim([-500 1500])
    ylim([0 0.5])
    title(num2str(wfq(freq),3))
    yl = ylim;
    hold on
    plot([-300 -300],[yl(1) yl(2)],'r--')
    plot([0 0],[yl(1) yl(2)],'k--')
    plot([700 700],[yl(1) yl(2)],'g--')
    hold off
    box off
end
subtitle('MRL')
%%
figure
hold on
for freq =  1:num_frequencies
    plot(-twin1:twin2-1,nanmean(all_mrls{freq}))
end
    plot([-300 -300],[yl(1) yl(2)],'r--')
    plot([0 0],[yl(1) yl(2)],'k--')
    plot([700 700],[yl(1) yl(2)],'g--')
    hold off
    box off
xlim([-500 1500])
%%
% save('CVTNEWLFPs_and_spikes')

%%
smval = 3;


degrees = [0:6:360]-180;
degrees = degrees*pi/180;

all_phases = cell(1,num_frequencies);
for freq = 1:num_frequencies
    all_phases{freq} = zeros(1,length(degrees));
end

all_mrls = cell(1,num_frequencies);
for unit = 1:size(phase_time,2)
    for freq = 1:num_frequencies
        phase = phase_time{freq,unit}(2,:);
        phase(isnan(phase)) = [];
         time = phase_time{freq,unit}(1,:);
         time(isnan(phase)) = []; %not sure why phase is NaN
         phase(isnan(phase)) = [];
         time(isnan(time)) = [];
%

        binned_phase = NaN(1,length(degrees));
        for bins = 2:length(degrees)
            these_dirs = (phase  < degrees(bins) & phase  >= degrees(bins-1));
            binned_phase(bins) = sum(these_dirs);
        end
        mrl = circ_r(phase);

        if mrl > 0.2 && length(phase) > 100
            plot(time,phase,'.')
            title(num2str(wfq(freq)))
            disp('now')
        end

        all_phases{freq} = all_phases{freq}+binned_phase;
        all_mrls{freq} = [all_mrls{freq} mrl];

    end
end
%%
figure
plot(wfq,cellfun(@median,all_mrls))
xlabel('Frequency (hz)')
ylabel('Median MRL')

%%
figure
for freq = 2:2:num_frequencies
    binned_phase = all_phases{freq};
    bf = [binned_phase(end-(3*smval):end) binned_phase binned_phase(1:3*smval)];%so don't get edege artifacts
    binned_phase = nandens(bf,smval,'gauss',1); %smooth to display firing rate by saccade direction
    binned_phase = binned_phase(3*smval+2:end-3*smval);%remove buffers

    subplot(4,4,freq/2)
    polar(degrees,binned_phase)
    title([num2str(wfq(freq),2) ' Hz']);

end

%

degrees = [0:1:360]-180;
degrees = degrees*pi/180;
for freq = 1:num_frequencies
    phase = phase_time{freq,unit}(2,:);
    time = phase_time{freq,unit}(1,:);
    time(isnan(phase)) = []; %not sure why phase is NaN
    phase(isnan(phase)) = [];
    time(isnan(time)) = [];

    binned_phase = NaN(1,length(degrees));
    for bins = 2:length(degrees)
        these_dirs = (phase  < degrees(bins) & phase  >= degrees(bins-1));
        binned_phase(bins) = sum(these_dirs);
    end
    mrl = circ_r(phase);

    bf = [binned_phase(end-(3*smval):end) binned_phase binned_phase(1:3*smval)];%so don't get edege artifacts
    binned_phase = nandens(bf,smval,'gauss',1); %smooth to display firing rate by saccade direction
    binned_phase = binned_phase(3*smval+2:end-3*smval);%remove buffers

    figure

    subplot(1,2,1)
    polar(degrees,binned_phase)
    title([num2str(wfq(freq)) ' Hz, MRL: ' num2str(mrl)]);

   subplot(1,2,2)
   plot(time,phase,'.')

end

%
clar
task = 'cvtnew';


fc = 6; % Cut off frequency
fs = 1000; % Sampling rate
[b,a] = butter(6,fc/(fs/2),'high'); % Butterworth filter of order 6

Cortex codes
trial_start_code = 15;
cross_on_code = 35;
fixation_code = 8;
dot_on_code = 25;
dot_clrchng_code = 27;
response_code = 4;
reward_code = 3;

min_blks = 1;

twin1 = 1000;
twin2 = 2000;

spike times vs spike phase
num_frequencies = 32;

all_LFPs = [];
filtered_all_LFPs = [];
for monkey = 2:-1:1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ---Read in Excel Sheet for Session data---%%%
    only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
        listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Vivian.mat'])
        
        predict_rt = 155;%155.85 ms prediction 5-percentile
        chamber_zero = [13.5 -11]; %AP ML
        
    elseif monkey ==2%strcmpi(monkey,'Tobii')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Tobii\';
        excel_file = [excel_dir 'Tobii_recordingnotes.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Recording Files\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\TO Figures\';
        
        predict_rt = 135;%ms prediction 5-percentile
        chamber_zero = [7.5 15]; %AP ML, his posertior hippocampus appears slightly shorter/more compressed than atlas
        
                listsq_read_excel(data_dir,excel_file);
        load([data_dir 'Across_Session_Unit_Data_Tobii.mat'])
        session_data(end) = [];%last file doesn't have strobe signal working so have no timing singnal :(
    end
    
    for sess = 1:length(session_data)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---import task and unit data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [task_file,~,~,multiunit,unit_names,unit_confidence,...
            sorting_quality,~,lfp_quality,~] = get_task_data(session_data{sess},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg','valid_trials','hdr');
        
        grap LPF channel info data
        num_trials = length(cfg.trl); %number of trials
        
        LFPchannels = NaN(1,4);
        try
            bad_channels = [];
            for channel = 1:4
                if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                    if  lfp_quality(channel) == 0; %if it is bad
                        bad_channels = [bad_channels channel];
                    else
                        temp = strfind(cfg.channel,['AD0' num2str(channel)]);
                        temp(cellfun(@isempty,temp)) = {0};
                        LFPchannels(channel) = find(cell2mat(temp));
                    end
                end
            end
            LFPchannels(bad_channels) = NaN;
        catch
            continue
        end
        
        if isempty(LFPchannels)
            continue
        end
        
        trial_LFP = cell(1,length(LFPchannels));
        filtered_trial_LFP = cell(1,length(LFPchannels));
        for chan = 1:length(LFPchannels);
            if ~isnan(LFPchannels(chan))
                trial_LFP{chan} = NaN(num_trials,twin1*2+twin2);
                filtered_trial_LFP{chan} = NaN(num_trials,twin1*2+twin2);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%---import & reformat data so that spikes are locked to events---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        disp([task_file(1:8) ' Aligning spike times to trial events'])
        
        num_trials = length(cfg.trl);
        
        for t = 1:num_trials
            if any(cfg.trl(t).allval == reward_code); %in which trial was successful
                
                trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == trial_start_code); %start of ITI
                trial_end = cfg.trl(t).alltim(end)-trial_start; %end of trial after reward
                doton = cfg.trl(t).alltim(cfg.trl(t).allval == dot_on_code)-trial_start; %when image turns on
                dotchange = cfg.trl(t).alltim(cfg.trl(t).allval == dot_clrchng_code)-trial_start; %when dot changes color
                monkey_response =  cfg.trl(t).alltim(cfg.trl(t).allval == response_code)-trial_start; %when release bar
                reward = cfg.trl(t).alltim(cfg.trl(t).allval == reward_code)-trial_start; %reward pulses
                
                                if dotchange-doton > 2000
                                    dotchange = doton+1999;
                                end
                trial_dur(t) = dotchange-doton;
                
                for chan = 1:length(LFPchannels)
                    if ~isnan(LFPchannels(chan))
                        LFPs = data(LFPchannels(chan)).values{t}; %spike trains for this trial
                        LFPs = LFPs(doton-twin1:dotchange+twin1);
                        trial_LFP{chan}(t,1:twin1+trial_dur(t)) =  LFPs(1:twin1+trial_dur(t));
                        filtered_trial_LFP{chan}(t,1:twin1+trial_dur(t)) =  filtfilt(b,a,LFPs(1:twin1+trial_dur(t)));
                    end
                end
                
                
                
                for chan = 1:length(LFPchannels)
                    if ~isnan(LFPchannels)
                        all_LFPs = [all_LFPs; nanmean(trial_LFP{chan})];
                        filtered_all_LFPs = [filtered_all_LFPs; nanmean(filtered_trial_LFP{chan})];
                    end
                end
            end
        end
    end
end
save('CVTNEWLFPs_filtered')
%
%%
figure
for freq = 2:2:num_frequencies
    subplot(4,4,freq/2)
    plot(tm,nanmean(all_LFP_phases{freq}))
    xlim([-500 twin2])
    title(num2str(wfq(freq),3))
    yl = ylim;
    hold on
    plot([-300 -300],[yl(1) yl(2)],'r--')
    plot([0 0],[yl(1) yl(2)],'k--')
    plot([700 700],[yl(1) yl(2)],'g--')
    hold off
    box off
end
subtitle('Mean Phase')
%%
figure
for freq =  2:2:num_frequencies
    subplot(4,4,freq/2)
    plot(-twin1:twin2-1,nanmean(all_mrls{freq}))
    xlim([-500 1500])
    ylim([0 0.5])
    title(num2str(wfq(freq),3))
    yl = ylim;
    hold on
    plot([-300 -300],[yl(1) yl(2)],'r--')
    plot([0 0],[yl(1) yl(2)],'k--')
    plot([700 700],[yl(1) yl(2)],'g--')
    hold off
    box off
end
subtitle('MRL')
%%
figure
hold on
for freq =  1:num_frequencies
    plot(-twin1:twin2-1,nanmean(all_mrls{freq}))
end
    plot([-300 -300],[yl(1) yl(2)],'r--')
    plot([0 0],[yl(1) yl(2)],'k--')
    plot([700 700],[yl(1) yl(2)],'g--')
    hold off
    box off
xlim([-500 1500])
%%
% save('CVTNEWLFPs_and_spikes')

%%
smval = 3;


degrees = [0:6:360]-180;
degrees = degrees*pi/180;

all_phases = cell(1,num_frequencies);
for freq = 1:num_frequencies
    all_phases{freq} = zeros(1,length(degrees));
end

all_mrls = cell(1,num_frequencies);
for unit = 1:size(phase_time,2)
    for freq = 1:num_frequencies
        phase = phase_time{freq,unit}(2,:);
        phase(isnan(phase)) = [];
         time = phase_time{freq,unit}(1,:);
         time(isnan(phase)) = []; %not sure why phase is NaN
         phase(isnan(phase)) = [];
         time(isnan(time)) = [];
%

        binned_phase = NaN(1,length(degrees));
        for bins = 2:length(degrees)
            these_dirs = (phase  < degrees(bins) & phase  >= degrees(bins-1));
            binned_phase(bins) = sum(these_dirs);
        end
        mrl = circ_r(phase);

        if mrl > 0.2 && length(phase) > 100
            plot(time,phase,'.')
            title(num2str(wfq(freq)))
            disp('now')
        end

        all_phases{freq} = all_phases{freq}+binned_phase;
        all_mrls{freq} = [all_mrls{freq} mrl];

    end
end
%%
figure
plot(wfq,cellfun(@median,all_mrls))
xlabel('Frequency (hz)')
ylabel('Median MRL')

%%
figure
for freq = 2:2:num_frequencies
    binned_phase = all_phases{freq};
    bf = [binned_phase(end-(3*smval):end) binned_phase binned_phase(1:3*smval)];%so don't get edege artifacts
    binned_phase = nandens(bf,smval,'gauss',1); %smooth to display firing rate by saccade direction
    binned_phase = binned_phase(3*smval+2:end-3*smval);%remove buffers

    subplot(4,4,freq/2)
    polar(degrees,binned_phase)
    title([num2str(wfq(freq),2) ' Hz']);

end

%

degrees = [0:1:360]-180;
degrees = degrees*pi/180;
for freq = 1:num_frequencies
    phase = phase_time{freq,unit}(2,:);
    time = phase_time{freq,unit}(1,:);
    time(isnan(phase)) = []; %not sure why phase is NaN
    phase(isnan(phase)) = [];
    time(isnan(time)) = [];

    binned_phase = NaN(1,length(degrees));
    for bins = 2:length(degrees)
        these_dirs = (phase  < degrees(bins) & phase  >= degrees(bins-1));
        binned_phase(bins) = sum(these_dirs);
    end
    mrl = circ_r(phase);

    bf = [binned_phase(end-(3*smval):end) binned_phase binned_phase(1:3*smval)];%so don't get edege artifacts
    binned_phase = nandens(bf,smval,'gauss',1); %smooth to display firing rate by saccade direction
    binned_phase = binned_phase(3*smval+2:end-3*smval);%remove buffers

    figure

    subplot(1,2,1)
    polar(degrees,binned_phase)
    title([num2str(wfq(freq)) ' Hz, MRL: ' num2str(mrl)]);

   subplot(1,2,2)
   plot(time,phase,'.')

end