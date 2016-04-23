m= 1/10*ones(1,4);
k(1) = 1/4; %4 Hz
k(2) = 1/16; %16 Hz
k(3) = 1/32; %30 Hz
k(4) =  1/64;%50 Hz
k = (2*pi./k).^2;
k = k.*m;
c = [0.25 1.5 0.5 0.5]; %larger more dampening, smaller less dampening
amplitude = [1 3 1 1];

data_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\';
figure_dir = 'C:\Users\seth.koenig\Documents\MATLAB\LFPAnalysis_anddecoding\';

load([data_dir 'BehavioralSaccadeRateData.mat'],'all_rates');
all_rates(isnan(all_rates)) = [];
all_rates = sort(all_rates);%'essentially cumulative distribution
num_rates = length(all_rates);

time = 0:0.001:30;% Define time vector
Ft = zeros(1,length(time)); %function over time

% %saccade type drive
% last_sac = 1;
% sacindex = [];
% while last_sac/1000 < max(time)-.8
%     rr = randi(num_rates);
%     last_sac = last_sac+all_rates(rr);
%     Ft(last_sac:last_sac+9) = 500; %arbitrary number
%     sacindex = [sacindex last_sac];
% end

%period drive e.g. model of rodents
drive_freq = 1000/10;
sacindex = [];
for t =1:length(time)/drive_freq-1
   Ft(drive_freq*t:drive_freq*t+9) = 500; 
   sacindex = [sacindex drive_freq*t];
end

Y = cell(1,4);
twin = 500;
avgSTA = zeros(1,2*twin+1);
titl = {'4','16','32','64','All'};

figure
for freq = 1:4
    disp(['Running for Freq# ' num2str(freq)])
    [~,Y] = ode45(@(t, y) HarmonicFcn(t,y,m(freq),k(freq),c(freq),Ft,time),time, [0, 1]);
    
    
    STA = zeros(1,2*twin+1);
    total_sacs = 0;
    
    for s = 1:length(sacindex)
        if sacindex(s) > twin && sacindex(s) < length(time)-twin
            STA = STA+Y(sacindex(s)-twin:sacindex(s)+twin,1)';
            total_sacs = total_sacs+1;
        end
    end
    
    STA = amplitude(freq)*STA/total_sacs;
    
    subplot(2,3,freq)
    hold on
    plot(-twin:twin,STA)
    plot([0 0],[min(STA) max(STA)],'k--')
    hold off
    grid on
    xlabel('Time')
    ylabel('"LFP"')
    title(titl{freq})
    
    avgSTA = avgSTA + STA;
end

subplot(2,3,5:6)
hold on
plot(-twin:twin,avgSTA)
plot([0 0],[min(avgSTA) max(avgSTA)],'k--')
hold off
grid on
xlabel('Time')
ylabel('"LFP"')
title(titl{end})
