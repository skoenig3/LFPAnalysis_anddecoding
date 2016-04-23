% Written to determine if LFPs are obviously locked to saccades
% Written November 21, 2014 by Seth  Konig

dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
fixation_durations = {};

cd(dir);
a = what;
m = a.mat;
set = 0;
average_LFP = cell(1,4); %set by electrode by sequence vs image
for file = 1:size(m,1);
    if ~isempty(strfind(m{file},'preprocessed'))
        fixationstats = [];
        load(m{file},'fixationstats','cfg','item_set','data','hdr');
        if ~isempty(strfind(item_set,'ckwenzr'))
           continue; 
        end
        lfpchans = find_desired_channels(hdr,data,'LFP');
        if ~isempty(fixationstats) %so is ListsQ
            set = set +1;
            for l = 1:4
%                 for ll = 1:2
                    average_LFP{set,l}= zeros(1,1000);
%                 end
            end
            [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
            [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
            num_sacs = 0;
            for t = 1:length(fixationstats);
                if any(cfg.trl(t).allval == 23) && itmlist(cfg.trl(t).cnd-1000) > sequence_items(end) %only want image trials
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
                    imgon =  cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start;%onlg get saccades while image is up
                    imgoff =  cfg.trl(t).alltim(cfg.trl(t).allval == 24)-trial_start;%
                    dataend =  cfg.trl(t).alltim(cfg.trl(t).allval == 20)-trial_start;%to make sure I have enough data for 500 ms after saccade
                    saccadetimes = fixationstats{t}.saccadetimes;
                    if ~isempty(saccadetimes)
                        tooearly = find(saccadetimes(1,:) <= imgon+500);%also want to ignore the 1st saccade since evoked by image
                        saccadetimes(:,tooearly) = [];
                        toolate = find(saccadetimes(1,:) >= imgoff);
                        saccadetimes(:,toolate) = [];
                        tootoolate =  find(saccadetimes(1,:) >= dataend-500);
                        saccadetimes(:,tootoolate) = [];
                        for s = 1:size(saccadetimes,2)
                            for l = 1:4
                                average_LFP{set,l} = average_LFP{set,l}+data(lfpchans(l)).values{t}(saccadetimes(1,s)-500:saccadetimes(1,s)+499);
                            end
                        end
                        num_sacs = num_sacs + size(saccadetimes,2);
                    end
                end
            end
            for l = 1:4
                average_LFP{set,l} = average_LFP{set,l}/num_sacs;
            end
        end
    end
end

figure, hold on
all_averages = zeros(1,1000);
for s = 1:6
    for l = 1:4
        plot(average_LFP{s,l})
        all_averages =all_averages+average_LFP{s,l};
    end
end

figure
plot(all_averages);
yl = ylim;
hold on
plot([501 501],[yl(1) yl(2)],'k')
xlabel('Time from Saccade (ms)')
ylabel('Normalized LFP Amplitude')
title('Average LFP averaged around the time of a Saccade')

%%
dir = 'C:\Users\seth.koenig\Documents\MATLAB\ListSQ\Recording Files\';
fixation_durations = {};

cd(dir);
a = what;
m = a.mat;
set = 0;
average_LFP = cell(1,4); %set by electrode by sequence vs image
for file = 1:size(m,1);
    if ~isempty(strfind(m{file},'preprocessed'))
        fixationstats = [];
        load(m{file},'fixationstats','cfg','item_set','data');
        lfpchans = length(data)-5:length(data)-2;
        if ~isempty(fixationstats) %so is ListsQ
            set = set +1;
            for l = 1:4
%                 for ll = 1:2
                    average_LFP{set,l}= zeros(1,1000);
%                 end
            end
            [itmlist,sequence_items,sequence_locations] = read_ListSQ_itm_and_cnd_files(item_set);
            [which_img,novel_vs_repeat] = get_image_numbers(cfg,itmlist,sequence_items,23);
            num_sacs = 0;
            for t = 1:length(fixationstats);
                if any(cfg.trl(t).allval == 23) && any(cfg.trl(t).allval == 3) % only want sequence trials
                    trial_start = cfg.trl(t).alltim(cfg.trl(t).allval == 15);
                    itemon =  cfg.trl(t).alltim(cfg.trl(t).allval == 23)-trial_start;%onlg get saccades while image is up
                    dataend =  cfg.trl(t).alltim(cfg.trl(t).allval == 20)-trial_start;%to make sure I have enough data for 500 ms after saccade
                    if length(dataend) > 1;
                        dataend = dataend(end);
                    end
                    saccadetimes = fixationstats{t}.saccadetimes;
                    if ~isempty(saccadetimes)
                        tooearly = find(saccadetimes(1,:) <= itemon+500);%also want to ignore the 1st saccade since evoked by image
                        saccadetimes(:,tooearly) = [];
                        tootoolate =  find(saccadetimes(1,:) >= dataend-500);
                        saccadetimes(:,tootoolate) = [];
                        for s = 1:size(saccadetimes,2)
                            for l = 1:4
                                average_LFP{set,l} = average_LFP{set,l}+data(lfpchans(l)).values{t}(saccadetimes(1,s)-500:saccadetimes(1,s)+499);
                            end
                        end
                        num_sacs = num_sacs + size(saccadetimes,2);
                    end
                end
            end
            for l = 1:4
                average_LFP{set,l} = average_LFP{set,l}/num_sacs;
            end
        end
    end
end
%%
figure, hold on
all_averages = zeros(1,1000);
for s = 1:6
    for l = 1:4
        plot(average_LFP{s,l})
        all_averages =all_averages+average_LFP{s,l};
    end
end

figure
plot(all_averages);
yl = ylim;
hold on
plot([501 501],[yl(1) yl(2)],'k')
xlabel('Time from Saccade (ms)')
ylabel('Normalized LFP Amplitude')
title('Average LFP averaged around the time of a Saccade')