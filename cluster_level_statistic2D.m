function [all_mc_matrix,pos_big_block_matrix,neg_big_block_matrix,observed_prctile] = cluster_level_statistic2D(observed_diff,shuffled_diff,min_cluster_area)
%written by Seth Konig March 10, 2017
%code implements method described by Maris and Oostenveld, 2007 to correct
%for multiple comparisons in 2D data. In particular time/frequency data comparing
%data across 2 difference conditions.

% Inputs
%   1) observed_diff: Observed 2D difference across conditions
%   2) shuffled_diff: Shuffled 2D differences
%   3) min_cluster_area: minimum cluster size to count as "significant"
%   cluster
%
% Outputs
%   1) all_mc_matrix: significant indeces after correcting for mulitple
%   comparison
%   2) pos_big_block_matrix: significant indeces before multiple comparions
%   corrections for positive differences across conditions
%   3) neg_big_block_matrix: significant indeces before multiple comparions
%   corrections for negative differences across conditions


alpha = 5;
if nargin < 3
    min_cluster_area =250;%area so e.g. Hz*ms
end
numshuffs = size(shuffled_diff,2);
num_time_pts = size(shuffled_diff{1},2);
num_frequncies = size(shuffled_diff{1},1);

%---Get Distrubution of Shuffled Values for Each Element in Time/Frequency---%
shuffled_distrubtion = NaN(num_frequncies,num_time_pts,numshuffs);
for shuff = 1:numshuffs
    shuffled_distrubtion(:,:,shuff) = shuffled_diff{shuff};
end

%---Get Observed Percentile---%
observed_prctile = NaN(num_frequncies,num_time_pts);
for freq = 1:num_frequncies
    for t = 1:num_time_pts
        shuffdist = shuffled_distrubtion(freq,t,:);
        observed_prctile(freq,t) = 100*sum(observed_diff(freq,t) > shuffdist(:))/numshuffs;
    end
end

%---Find Significant Clusters---%

%find "positive" differences
pos_potential_sig_times = (observed_prctile > 97.5);
pos_CC = bwconncomp(pos_potential_sig_times);
pos_block_size = cellfun(@numel,pos_CC.PixelIdxList);
pos_big_blocks = find(pos_block_size > min_cluster_area);
pos_big_block_matrix = zeros(size(pos_potential_sig_times));
pos_block_difference = zeros(1,length(pos_big_blocks));
for b = 1:length(pos_big_blocks)
    pos_big_block_matrix(pos_CC.PixelIdxList{pos_big_blocks(b)}) = 1;
    pos_block_difference(b) = sum(sum(observed_diff(pos_CC.PixelIdxList{pos_big_blocks(b)})));
end

%find "negative" differences
neg_potential_sig_times = (observed_prctile < 2.5);
neg_CC = bwconncomp(neg_potential_sig_times);
neg_block_size = cellfun(@numel,neg_CC.PixelIdxList);
neg_big_blocks = find(neg_block_size > min_cluster_area);
neg_big_block_matrix = zeros(size(neg_potential_sig_times));
neg_block_difference = zeros(1,length(neg_big_blocks));
for b = 1:length(neg_big_blocks)
    neg_big_block_matrix(neg_CC.PixelIdxList{neg_big_blocks(b)}) = 1;
    neg_block_difference(b) = sum(sum(observed_diff(neg_CC.PixelIdxList{neg_big_blocks(b)})));
end

%---Find Biggest Cluster---%
%abs/magnitude of difference
if isempty(pos_block_difference) && ~isempty(neg_block_difference)
    [~,biggest_cluster] = min(neg_block_difference);
    biggest_cluster_ind = neg_CC.PixelIdxList{neg_big_blocks(biggest_cluster)};
elseif ~isempty(pos_block_difference) && isempty(neg_block_difference)
    [~,biggest_cluster] = max(pos_block_difference);
    biggest_cluster_ind = pos_CC.PixelIdxList{pos_big_blocks(biggest_cluster)};
elseif max(abs(neg_block_difference)) > max(pos_block_difference)
    [~,biggest_cluster] = min(neg_block_difference);
    biggest_cluster_ind = neg_CC.PixelIdxList{neg_big_blocks(biggest_cluster)};
else
    [~,biggest_cluster] = max(pos_block_difference);
    biggest_cluster_ind = pos_CC.PixelIdxList{pos_big_blocks(biggest_cluster)};
end

%define new shuffled difference distrubtion as distribution of biggest cluster
new_shuffled_dist = NaN(1,numshuffs);
for shuff = 1:numshuffs
    new_shuffled_dist(shuff) = sum(shuffled_diff{shuff}(biggest_cluster_ind));
end


%---Determine Significant Clusters after Correcting for Multiple Comparisons---%
neg_block_difference_percentile = NaN(1,length(neg_block_difference));
pos_block_difference_percentile = NaN(1,length(pos_block_difference));
all_mc_matrix = zeros(size(pos_potential_sig_times));
for n = 1:length(neg_block_difference);
    neg_block_difference_percentile(n) = 100*sum(neg_block_difference(n) > new_shuffled_dist)/numshuffs;
    if neg_block_difference_percentile(n) < 2.5
        all_mc_matrix(neg_CC.PixelIdxList{neg_big_blocks(n)}) = 1;
    end
end
for p = 1:length(pos_block_difference)
    pos_block_difference_percentile(p) = 100*sum(pos_block_difference(p) > new_shuffled_dist)/numshuffs;
    if  pos_block_difference_percentile(p) > 97.5
        all_mc_matrix(pos_CC.PixelIdxList{pos_big_blocks(p)}) = 1;
    end
end

end