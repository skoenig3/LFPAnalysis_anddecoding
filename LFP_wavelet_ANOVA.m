dbqfunction contrasts = LFP_wavelet_ANOVA(parsed_LFP_data,factors,factortype)
% adapted for LFP analysis by Seth Konig August, 2014
% wfANOVAdemo.m
%
% sample code that recreates the primary analysis figure of
% "Statistically-significant contrasts between EMG waveforms revealed using wavelet-based functional ANOVA"
% McKay, Welch, Vidakovic, and Ting
% in review, Journal of Neurophysiology
%
% usage: type "wfANOVAdemo"
%
% modifications or playing around with the statistics or wavelet transform
% (in the "wfANOVA" subfunction) should not generally affect the plotting
% of contrasts, but may affect the big plot that gets generated - so both
% are provided.
%
% please contact me (j.lucas.mckay@emory.edu) if you have any questions or
% difficulties running the code!
%
% J. Lucas McKay, 19 August 2012

% load data for tibialis anterior
% load wfANOVAdata
% data = TA;
% factors = [V A S];

% the data are as follows:
%  Tibialis anterior EMG, replicates x time samples:
%   TA        439x512            1798144  double
%  Time vector:
%   time        1x512               4096  double
%  Factor vectors, coding for peak velocity, peak acceleration, and
%  subject:
%   V         439x1                 3512  double
%   A         439x1                 3512  double
%   S         439x1                 3512  double
if nargin < 3
    factortype = [];
end

time = 1:size(parsed_LFP_data,2);%in ms since sampled at 1000 Hz

% perform post-hoc tests on factors V and A:
performposthoc = 1;

% calculate contrasts using wavelet-based functional ANOVA.
[contrasts] = wfANOVA(parsed_LFP_data,factors,performposthoc);

% plot a large figure to match the primary analysis figure in the paper.
plotbigfigure(parsed_LFP_data,time,factors,contrasts,factortype)

end

function varargout = wfANOVA(data,factors,performposthoc)

% wavelet transform the data
[wavedata, waveparams] = wavelettransform(data);

% to test for perfect reconstruction, you can use the following:
% recondata = inversewavelettransform(wavedata,waveparams);
alpha = 0.05;
display = 'off'; % display ANOVA tables? on/off
% Set up number of tests and length of records.
[ntrls,npnts] = size(data);
% change the grouping variable to a cell array.
group = {};
for j = 1:size(factors,2)
    group{1,j} = factors(:,j);
end
% perform pointwise ANOVA
t_p = [];
stats = {};
han = waitbar(0,'performing pointwise ANOVA','toolbar','none');
for j = 1:npnts
    [t_p(:,j),~,stats{j}] = anova1(wavedata(:,j),factors,display,alpha);
    waitbar(j/npnts,han);
end
close(han)
contrasts = {};
for j = 1:length(performposthoc)
    % perform post-hoc tests for each flagged factor.
    if performposthoc(j)
        tempwavecontrast = posthocsub(t_p(j,:),alpha,stats,j);
        temptimecontrast = inversewavelettransform(tempwavecontrast,waveparams);
        contrasts{j} = temptimecontrast;
    end
end
varargout = contrasts(logical(performposthoc));
end

function contrasts = posthocsub(pvals,alpha,stats,dim)
% type of critical value used for post hoc tests
ctype = 'scheffe';
% don't display detailed results
display = 'off';
% calculate p-values for post-hoc tests
posthocalpha = alpha/sum(pvals<alpha);
% figure out the maximum levels for post-hoc tests from the stats
% structure
maxcontrast = length(stats{1}.gnames)-1;
% assemble the contrast structure
contrasts = zeros(maxcontrast,length(pvals));
% loop through and do the comparison
for j = 1:length(stats)
    if pvals(j)<alpha
        clear temp
        [temp.contrasts, temp.means] = multcompare(stats{j},'dimension',dim,'display',display,'ctype',ctype,'alpha',posthocalpha);
        for k = 1:size(temp.contrasts,1)
            if (sign(temp.contrasts(k,3)) == sign(temp.contrasts(k,5)))
                % Include the contrast in the contrast waveform. Note the
                % negative sign is because Matlab outputs level 1 - level
                % 2, and we want level 2 - level 1.
                contrasts(k,j) = -temp.contrasts(k,4);
            else
                contrasts(k,j) = 0;
            end
        end
    end
end
end

function [wavedata, waveparams] = wavelettransform(data)
% Get the length of the data records
[ntrls,npnts] = size(data);
% Wavelet to use: 3-rd order Coiflet.
wavestr = 'coif3';
% Extension mode: periodic, so that what is returned is the same length as
% the data.
wavemode = 'per';
dwtmode(wavemode);
% Find the maximum possible level of the wavelet decomposition:
lev = wmaxlev([1 npnts],wavestr);
% Loop through and transform
wavedata = zeros(size(data));
for i = 1:ntrls
    x = data(i,:);
    [wx,wL] = wavedec(x,lev,wavestr);
    wavedata(i,:) = wx;
end
waveparams.wL = wL;
waveparams.wavestr = wavestr;
waveparams.wavemode = wavemode;
end

function data = inversewavelettransform(wavedata,waveparams)
dwtmode(waveparams.wavemode);
% Get the length of the data records
[ntrls,npnts] = size(wavedata);
data = zeros(size(wavedata));
for i = 1:ntrls
    data(i,:) = waverec(wavedata(i,:),waveparams.wL,waveparams.wavestr);
end
end

function h = plotbigfigure(parsed_LFP_data,time,factors,contrasts,factortype)
if max(factors) == 2;
    figure
    % plot the differences in LFPs between Sequences:
    subplot(3,1,1)
    plot(mean(parsed_LFP_data(factors == 1,:)));
    grid on
    title(['Mean LFP for Factor 1 n = ',num2str(sum(factors == 1))])
    if max(time) == 512
        xlim([12 512])
    else
        xlim([0 max(time)])
    end
    xlabel('Time (ms)')
    ylabel('Normalized LFP (uV)')
    subplot(3,1,2)
    plot(mean(parsed_LFP_data(factors == 2,:)));
    grid on
    if max(time) == 512
        xlim([12 512])
    else
        xlim([0 max(time)])
    end
    xlabel('Time (ms)')
    ylabel('Normalized LFP (uV)')
    title(['Mean LFP for Factor 2 n = ',num2str(sum(factors == 1))])
    subplot(3,1,3);
    hold on
    plot(mean(parsed_LFP_data(factors == 2,:))-mean(parsed_LFP_data(factors == 1,:)));
    plot(contrasts,'r')
    hold off
    grid on
    if max(time) == 512
        xlim([12 512])
    else
        xlim([0 max(time)])
    end
    xlabel('Time (ms)')
    ylabel('Normalized LFP (uV)')
    legend('Mean Difference','wfANOVA')
elseif max(factors == 4);
    figure
    % plot the differences in LFPs between recording channels or items in
    % sequence
    contrast_count = 1;
    [i,j] = find(tril(ones(4,4)));
    for channels = 1:length(i);
        subplot(4,4,sub2ind([4,4],i(channels),j(channels)))
        if i(channels) == j(channels) %just plot mean LFP
            plot(mean(parsed_LFP_data(factors == i(channels),:)));
            grid on
            title(['Mean LFP ' factortype ' ' num2str(i(channels))])
        else
            hold on
            plot(mean(parsed_LFP_data(factors == i(channels),:)-parsed_LFP_data(factors == j(channels),:)));
            plot(contrasts(contrast_count,:),'r')
            hold off
            grid on
            title(['LFP differences on ' factortype ' ' num2str(j(channels)) ' and ' num2str(i(channels))])
            contrast_count = contrast_count + 1;
        end
        xlabel('Time (ms)')
        ylabel('LFP (uV)')
        xlim([12 500])
    end
end
end