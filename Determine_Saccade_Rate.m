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


saccade_rate = cell(1,length(listsq_files));


for file = 1:length(listsq_files)
    disp(['File #' num2str(file) ': ' listsq_files{file}(1:end-11) ])
    
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
    
    inter_saccade_interval = NaN(num_trials,75);
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
                fixationtimes(:,toolate) = NaN;
                
                tooearly = find(saccadetimes(1,:) <= imgon-trial_start+250);
                saccadetimes(:,tooearly) = NaN;
                tooearly = find(fixationtimes(1,:) <= imgon-trial_start+250);
                fixationtimes(:,tooearly) = NaN;
                
                %time outside image is ignored so not 1:1 saccaddes to
                %fixations
                fixdur =  NaN(1,size(saccadetimes,2));
                for sac = 1:size(saccadetimes,2);
                    sacstart = saccadetimes(1,sac);
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
                
                %this may not be the best way but it is a way. Might be too
                %conservative
                for s = 1:length(goodsacs)-1
                    sac = goodsacs(s);
                    if sac+1 > size(fixationtimes,2)
                        continue
                    end
                    if fixationtimes(1,sac+sacindex)-saccadetimes(2,sac) == 1%there is a fixation on the image after the saccade
                        if saccadetimes(1,sac+1)-fixationtimes(2,sac+sacindex) == 1 %there is another saccade on the image after this
                            inter_saccade_interval(t,sac) = saccadetimes(1,sac+1)-saccadetimes(1,sac);
                            %start of saccade to start of next saccade;
                        end
                    end
                end
            end
        end
    end
    
    inter_saccade_interval = inter_saccade_interval(novel_vs_repeat == 1,:);
    inter_saccade_interval(isnan(inter_saccade_interval)) = [];
    saccade_rate{file} = inter_saccade_interval;
end

mean_rates = NaN(1,length(listsq_files));
median_rates = NaN(1,length(listsq_files));
std_rates = NaN(1,length(listsq_files));
all_rates = NaN(1,length(listsq_files)); 

for file = 1:length(listsq_files)
    mean_rates(file) = mean(saccade_rate{file});
    std_rates(file) = std(saccade_rate{file});
    median_rates(file) = median(saccade_rate{file});
    all_rates = [all_rates saccade_rate{file}];
end

figure
hist(1000./all_rates,100)
xlabel('Saccade Rate')
ylabel('Count')
title('Saccade Rate during Novel Images')
save_and_close_fig(figure_dir,'Histogram_of_SaccadeRates')
save('BehavioralSaccadeRateData')
