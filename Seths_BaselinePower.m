function BaseLinePower = Seths_BaselinePower(freq,time_interval)
%written by Seth Konig 4/28/16 
%code from fieldtrip just cut what I needed to can use any desired period
%    
% Inputs:
%   a) TFdata: power spectrum
%   b) time_interval of interest
%
% Outputs:
%   a) TFbl: baseline power by frequency (row channel, col is frequency)


TFdata = freq.powspctrm;
timeVec = freq.time;

tidx = find(timeVec >= time_interval(1) & timeVec <= time_interval(2));

TFbl = NaN(size(TFdata,1),size(TFdata,2));
for k=1:size(TFdata,2) % loop frequencies
    for l=1:size(TFdata,1) % loop channels
        TFbl(l,k) = nanmean(TFdata(l,k,tidx),3); %compute average baseline power
        if TFbl(l,k) == 0,
            error('Average baseline power is zero');
        end
    end
end
BaseLinePower = TFbl;
end