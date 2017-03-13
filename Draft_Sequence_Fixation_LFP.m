clar

task = 'ListSQ';
twin = 500;% how much time to take before and after saccade.
% Delay in HPC to visual stimuli around 150-200 ms so probably want at least 2x this
fixwin = 5;%size of fixation window on each crosshair
item_event_codes = [23 24 25 26 27 28 29 30]; %odd indeces item on, even indeces item off
trial_start_code = 15;
smval = 50;%gaussian 1/2 width for smoothing
numshuffs = 1000; %recommend this is between 100 & 1000, for bootstrapping to
% dtermine significance
info_type = 'temporal';
Fs = 1000;
min_blks = 2;


Fixation_trig_avg = [];
all_reaction_times = [];
all_orders = [];
all_monkeys = [];
for monkey = 1:2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %---Read in Excel Sheet for Session data---%%%
    %only need to run when somethings changed or sessions have been added
    if monkey == 1%strcmpi(monkey,'Vivian')
        excel_dir = '\\towerexablox.wanprc.org\Buffalo\eblab\PLX files\Vivian\';
        excel_file = [excel_dir 'Vivian_Recording_Notes-ListSQ.xlsx']; %recording notes
        data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted\';
        figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\PW Resorted Figures\';
        
%         listsq_read_excel(data_dir,excel_file);
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
        task = 'ListSQ';
        [task_file,item_file,cnd_file,multiunit,unit_names,unit_confidence,...
            sorting_quality,waveform_count,lfp_quality,comments] = get_task_data(session_data{sess},task);
        if isempty(task_file)
            warning('No file could be found for specificed task. Exiting function...')
            continue
        end
        
        %load trial data
        load([data_dir task_file(1:end-11) '-preprocessed.mat'],'data','cfg',...
            'hdr','fixationstats');
        num_trials = length(cfg.trl);
        disp(['Session #' num2str(sess)])
        
        LFPchannels = find_desired_channels(cfg,'LFP');
        %remove bad LFP channels
        bad_channels = [];
        for channel = 1:4
            if cell2mat(strfind(hdr.label,['AD0' num2str(channel)])) %make sure have recorded channel
                if  lfp_quality(channel) == 0; %if it is bad
                    bad_channels = [bad_channels channel];
                end
            end
        end
        try
            LFPchannels(bad_channels) = [];
        catch
            continue
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Get successful trials Information by Task---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %get important task specific information
        [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_file,cnd_file);
        overlap = find((sequence_locations{1}(1,:) ==  sequence_locations{2}(1,:)) & ...
            (sequence_locations{1}(2,:) ==  sequence_locations{2}(2,:)));
        
        
        %preallocate space and parallel structure of cfg
        successful_sequence_trials = NaN(1,length(cfg.trl));
        which_sequence = NaN(1,length(cfg.trl));
        for t = 1:num_trials;
            if sum(cfg.trl(t).allval == 3) >= 5; %in which sequence trials were rewarded
                which_sequence(t) = find(sequence_items == itmlist(cfg.trl(t).cnd-1000));
                successful_sequence_trials(t) = t;
            end
        end
        successful_sequence_trials = laundry(successful_sequence_trials);
        which_sequence = laundry(which_sequence);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---process eye data locked to trial events---%%%
        num_trials = length(successful_sequence_trials);
        fixstats = fixationstats;
        fixationstats = fixationstats(successful_sequence_trials);
        cfg.trl = cfg.trl(successful_sequence_trials);
        
        fixation_start_time = NaN(length(fixationstats),4);%when did fixation on item start
        
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
                event_codes,event_times,1);
            
            fixation_numbers = trialdata.fixationnums; %fixation number for each item
            fixationtimes = fixationstats{trial}.fixationtimes;
            for item = 1:4
                if ~isnan(fixation_numbers(item))
                    fixation_start_time(trial,item) = fixationtimes(1,fixation_numbers(item));
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%---Calculate Firing Rate Locked to Eye Data---%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %also locked to "trial start" item 1 on or saccade to item 1
        trial_LFPs = NaN(4*num_trials*length(LFPchannels),2*twin+1);
        reaction_time = NaN(1,4*num_trials);
        order = NaN(1,4*num_trials);
        which_monkey = NaN(1,4*num_trials);
        ind = 1;
        for trial = 1:num_trials
            for c = 1:4;
                fixt = fixation_start_time(trial,c);
                
                item_on = cfg.trl(trial).alltim(cfg.trl(trial).allval == item_event_codes(2*c-1)); %when item turned on
                item_on = item_on-cfg.trl(trial).alltim(1);
                
                if isnan(fixt)
                    continue
                elseif fixt < twin %starts before trial does
                    continue
                end
                
                
                for chan = 1:length(LFPchannels)
                    trial_LFPs(ind,:) = data(LFPchannels(chan)).values{successful_sequence_trials(trial)}(fixt-twin:fixt+twin);
                    reaction_time(ind) = fixt-item_on;
                    order(ind) = c;
                    which_monkey(ind) = monkey;
                    ind = ind+1;
                end
            end
        end
        trial_LFPs = laundry(trial_LFPs);
        reaction_time = laundry(reaction_time);
        order = laundry(order);
        which_monkey = laundry(which_monkey);
        
        Fixation_trig_avg = [Fixation_trig_avg; trial_LFPs];
        all_reaction_times = [all_reaction_times reaction_time];
        all_orders = [all_orders order];
        all_monkeys = [all_monkeys which_monkey];       
    end
end
%%
truly_predcit = Fixation_trig_avg(all_reaction_times < 0  & all_reaction_times > -500 & all_orders > 1,:);
predict = Fixation_trig_avg(all_reaction_times < 135  & all_reaction_times > 0 & all_orders > 1,:);
fast_react = Fixation_trig_avg(all_reaction_times < 200  & all_reaction_times > 155 & all_orders > 1,:);
modest_react = Fixation_trig_avg(all_reaction_times < 275  & all_reaction_times > 225 & all_orders > 1,:);
slow_react = Fixation_trig_avg(all_reaction_times < 400  & all_reaction_times > 300 & all_orders > 1,:);
%%
figure
hold on
plot(tm,mean(truly_predcit));
plot(tm,mean(predict));
plot(tm,mean(fast_react));
plot(tm,mean(modest_react));
plot(tm,mean(slow_react));
hold off
xlabel('Time From Saccde Start (ms)')
ylabel('LFP (uV)')
xlim([-200 500])

legend('Truly Predictive','Predictive','Fast Reactive','Modest Reactive','Slow Reative')

%%

[~,~,wfq,truly_predcit_meanpower,truly_predcit_meanphase,~,truly_predcit_phasevar] = waveletanalysis(truly_predcit(1:10:end,:));
[~,~,wfq,predict_meanpower,predict_meanphase,~,predict_phasevar] = waveletanalysis(predict(1:14:end,:));
%%
[~,~,wfq,fast_react_meanpower,fast_react_meanphase,~,fast_react_phasevar] = waveletanalysis(fast_react(1:100:end,:));
[~,~,wfq,modest_react_meanpower,modest_react_meanphase,~,modest_react_phasevar] = waveletanalysis(modest_react(1:60:end,:));
[~,~,wfq,slow_react_meanpower,slow_react_meanphase,~,slow_react_phasevar] = waveletanalysis(slow_react(1:8:end,:));
%%

tm = -500:500;
figure
subplot(2,2,1)
low = find(wfq < 30);
imagesc(tm,wfq(low),slow_react_meanpower(low,:))
xlabel('Time from  Fixation start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')

subplot(2,2,4)
mid = find(wfq >= 30 & wfq < 80);
imagesc(tm,wfq(mid),slow_react_meanpower(mid,:))
xlabel('Time from  Fixation start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')


subplot(2,2,2)
high = find(wfq > 70);
imagesc(tm,wfq(high),slow_react_meanpower(high,:))
xlabel('Time from  Fixation start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')
%%
imagesc(tm,wfq,truly_predcit_meanphase-predict_meanphase)
xlabel('Time from  Fixation start (ms)')
ylabel('Frequency (Hz)')
axis xy
colormap('jet')