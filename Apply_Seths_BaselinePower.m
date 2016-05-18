function freq = Apply_Seths_BaselinePower(freq,baseline_powerspectrum,baselinetype)
%written by Seth Konig 4/28/16 
%code from fieldtrip just cut what I needed to can use any desired period
%    
% Inputs:
%   a) freq: frequency structure with powerspectrum
%   b) time_interval of interest
%   c) baselinetypes: 'absolute','relchange','relative'
%
% Outputs:
%   a) freq: normalized frquency

TFdata = freq.powspctrm;
TFbl = baseline_powerspectrum;

if strcmpi(baselinetype,'relative')
    for k=1:size(TFdata,2) % loop frequencies
        for l=1:size(TFdata,1) % loop channels
            TFdata(l,k,:) = TFdata(l,k,:) / TFbl(l,k);     % compute relative change (i.e. ratio)
        end
    end
    
elseif strcmpi(baselinetype,'absolute')
    for k=1:size(TFdata,2) % loop frequencies
        for l=1:size(TFdata,1) % loop channels
            TFdata(l,k,:) = TFdata(l,k,:) - TFbl(l,k);        % subtract baseline power
        end
    end
    
elseif strcmpi(baselinetype,'relchange')
    for k=1:size(TFdata,2) % loop frequencies
        for l=1:size(TFdata,1) % loop channels
            TFdata(l,k,:) = ((TFdata(l,k,:) - TFbl(l,k)) / TFbl(l,k)); % compute relative change
        end
    end
end

freq.powspctrm = TFdata;