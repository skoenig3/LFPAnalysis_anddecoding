clar
data_dir = '\\BuffaloTower.wanprc.org\Buffalo\Miriam\SFT_preprocessedData\eyeDataFromNs6andConcatedWithLFP\';


fs = 1000;
[bhigh,ahigh] = butter(6,65/(fs/2),'high');


min_fix_dur = 150;
min_sac_amp = 2;

data_files = {'wr130212sft-cut_Images','wr130213sft-cut_Images','wr130214sft-cut_Images',...
    'wr130215sft-cut_Images','wr130218sft-cut_Images','wr130219sft-cut_Images','wr130220sft-cut_Images',...
    'wr130227sft-cut_Images','wr130301sft-cut_Images'};


twin = 500;
STA = cell(length(data_files),13);
for chan = 1:13
    for sess = 1:length(data_files)
        STA{sess,chan} = NaN(3500,twin*2+1);
    end
end

for sess = 1:length(data_files)
    disp(['Loading ' data_files{sess}])
    load([data_dir  data_files{sess}])
  
    
    cfg.channel = tl.timelock.label;
    [eye_channels] = find_desired_channels(cfg,'Eye');
    [lfp_channels] = find_desired_channels(cfg,'LFPs');

    
    channeldata = tl.timelock.trial;
    num_trials = size(channeldata,1);
    
    sacind = ones(1,length(lfp_channels));
    for t = 1:num_trials
        x = squeeze(channeldata(t,eye_channels(1),:));
        y = squeeze(channeldata(t,eye_channels(2),:));
        nan_ind = find(isnan(x),1,'first');
        x = x(1:nan_ind-1);
        y = y(1:nan_ind-1);
        xy = [x';y'];
        len = length(x);
        if len < 2500
            continue
        end
        
        [fixationtimes,saccadetimes] = adaptive_VT_Recording_Data(xy);

        for chan = 1:length(lfp_channels)
            LFPdata = squeeze(channeldata(t,lfp_channels(chan),:))';
            %LFPdata = filtfilt(bhigh,ahigh,LFPdata(1:len));
            for s = 1:size(saccadetimes,2)
                sact = saccadetimes(1,s);
                if sact > twin && sact < len-twin
                    sacend =  saccadetimes(2,s);
                    
                    post_fix = find(fixationtimes(1,:)-1 ==sacend);%next fixation should start immediately after
                    if isempty(post_fix) %trial ended or eye outside of image
                        continue %try next one
                    end
                    sacamp = sqrt(sum((xy(:,sacend)-xy(:,sact)).^2));%saccade amplitude
                    fix_dur = fixationtimes(2,post_fix)-fixationtimes(1,post_fix)+1;%this fixation duration
                    
                    if sacamp >= min_sac_amp && fix_dur >= min_fix_dur %next fixation has to be long enough & Fixation large enough                       
                        STA{sess,chan}(sacind(chan),:) = LFPdata(sact-twin:sact+twin);
                        sacind(chan) = sacind(chan)+1;
                    end
                end
            end
        end
    end
end
%%
average_STAs = NaN(length(data_files)*12,2*twin+1);
ind = 1;
for sess = 1:length(data_files)
   for chan = 1:12
       average_STAs(ind,:) = nanmean(STA{sess,chan});
       ind = ind+1;
   end
end

%%
ref_average_STAs = NaN(length(data_files),2*twin+1);
ind = 1;
for sess = 1:length(data_files)
   for chan = 13
       ref_average_STAs(ind,:) = nanmean(STA{sess,chan});
       ind = ind+1;
   end
end
%%
figure
subplot(1,2,1)
hold on
plot(-twin:twin,average_STAs');
plot([-twin twin],[0 0],'k')
yls = ylim;
yld = max(abs(yls));
plot([0 0],[-yld yld],'k--')
hold off
xlabel('Time From Saccade Onset (ms)')
ylabel('LFP (uV)')
xlim([-250 500])
title('Individual STAs')

subplot(1,2,2)
hold on
plot(-twin:twin,nanmean(average_STAs));
plot([-twin twin],[0 0],'k')
yls = ylim;
yld = max(abs(yls));
plot([0 0],[-yld yld],'k--')
hold off
xlabel('Time From Saccade Onset (ms)')
ylabel('LFP (uV)')
xlim([-250 500])
title('Average STAs')
