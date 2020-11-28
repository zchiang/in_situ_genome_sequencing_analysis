% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Repetitive element analyses for Fig. 2

%% Load full table

full_table = readtable('../data/Table_S1_pgp1_data_table.csv');

%% Get all unique repetitive elements

count_threshold = 50;

% find unique classes

[uniq_classes ia ic] = unique(full_table.rep_class);
class_counts = accumarray(ic,1);
filt_uniq_classes = uniq_classes(class_counts>count_threshold & ~contains(uniq_classes,"N/A") & ~contains(uniq_classes,"?") & string(uniq_classes)~='');
filt_class_counts = class_counts(class_counts>count_threshold & ~contains(uniq_classes,"N/A") & ~contains(uniq_classes,"?") & string(uniq_classes)~='');

% find unique families

[uniq_fams ia ic] = unique(full_table.rep_family);
fam_counts = accumarray(ic,1);
filt_uniq_fams = uniq_fams(fam_counts>count_threshold & ~contains(uniq_fams,"N/A") & ~contains(uniq_fams,"?") & string(uniq_fams)~='');
filt_fam_counts = fam_counts(fam_counts>count_threshold & ~contains(uniq_fams,"N/A") & ~contains(uniq_fams,"?") & string(uniq_fams)~='');

% find unique names

[uniq_names ia ic] = unique(full_table.rep_name);
name_counts = accumarray(ic,1);
filt_uniq_names = uniq_names(name_counts>count_threshold & ~contains(uniq_names,"N/A") & ~contains(uniq_names,"?") & string(uniq_names)~='');
filt_name_counts = name_counts(name_counts>count_threshold & ~contains(uniq_names,"N/A") & ~contains(uniq_names,"?") & string(uniq_names)~='');

% combine all annotation levels

all_uniq_elements = [filt_uniq_classes; filt_uniq_fams; filt_uniq_names];
all_counts = [filt_class_counts; filt_fam_counts; filt_name_counts];
all_levels = [repmat(2,size(filt_class_counts)); repmat(3,size(filt_fam_counts)); repmat(1,size(filt_name_counts))];

%% Calculate enrichment z-scores for each RE - bin combo

% set parameters

bin_size = 0.1;
step_size = 0.01;

range = 0:step_size:max_dist;
num_perms = 500;
num_bins = length(range);
table_size = size(full_table,1);

% create # bins by # RE matrix for z-scores

bin_zscore = zeros(num_bins,size(all_uniq_elements,1));
    
% loop through all unique names

for re_index=1:size(all_uniq_elements,1)
    
    % count occurences of repetitive element
    
    re = all_uniq_elements(re_index)
    re_freq = all_counts(re_index);

    % permute repetitive element labels
    
    background = zeros(num_perms,re_freq);
    for j=1:num_perms
        background(j,:) = full_table.norm_r_2D(randsample(table_size,re_freq));
    end

    % loop through bins
    
    bin_index = 1;
    for j=range
        
        % define bin boundaries
        
        start_dist = max(j-bin_size/2,0);
        end_dist = min(j+bin_size/2,1);

        % calcualte z-score
        
        exp = sum(start_dist <= background & background < end_dist,2);
        obs = sum(start_dist <= full_table.norm_r_2D & full_table.norm_r_2D < end_dist & string(full_table{:,22+all_levels(re_index)}) == re);

        bin_zscore(bin_index,re_index) = (obs-mean(exp))./std(exp);
        bin_index = bin_index + 1;
        
    end
    
end

%% Calculate variability

zscores_var = std(bin_zscore).^2;
[sort_var sort_order] = sort(zscores_var,'descend');
sort_elements = all_uniq_elements(sort_order);
sort_zscore = bin_zscore(:,sort_order);

% remove duplicate Simple_repeat and Satellite elements

sort_var(:,[4,10]) = [];
sort_elements([4,10]) = [];
sort_zscore(:,[4,10]) = [];

%% Plot REs by variability (Fig. 2F)

thresh = 1.75;
num_labels = sum(sort_var>thresh);

figure;
scatter(1:size(sort_var,2),sort_var,25,'k','.'); hold on;
scatter(1:num_labels,sort_var(1:num_labels),50,'r','.'); hold on;

text([1:num_labels]+2,sort_var(1:num_labels),strrep(sort_elements(1:num_labels),'_','\_'))
xlabel('Repetitive elements','FontSize',12); ylabel('Radial bias','FontSize',12)
x = linspace(-5,300); y = x./x*thresh;
plot(x,y,'--','Color',[0.7 0.7 0.7]); hold on;
xlim([0 size(sort_var,2)]); yticks([0:3:12])

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 3 3];  

%% Select most variable elements for heatmap and sort by max bin

most_var = sort_zscore(:,1:num_labels);
most_var_elements = sort_elements(1:num_labels);

[max_val max_bin] = max(most_var,[],1);
[sort_bin bin_order] = sort(max_bin);

%% Plot RE enrichment heatmap (Fig. 2G)

figure; imagesc(log10(normcdf(most_var(:,bin_order))).*~sign+(-1*(log10(normcdf(-1*most_var(:,bin_order))).*sign)))

%colormap(redblue); caxis([-5 5]);
colormap(colorsJDB(0,0,'solar_extra')); caxis([-5 5]);
h = colorbar; xlabel(h, 'signed log10(p-value)')

yticks([1:25:101]); yticklabels({string(0:0.25:1)})
xticks([1:num_labels])
xticklabels({strrep(string(most_var_elements(bin_order)),"_","\_")})
xtickangle(90)
set(gca,'fontsize',10);
ylabel({'Rel. radial distance','from center'})

ax = gca; ax.YDir = 'normal'

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 5 3];  