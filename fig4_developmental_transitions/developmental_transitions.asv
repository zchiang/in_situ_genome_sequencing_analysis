% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Spatial separation, Rabl, and GC analyses for Fig. 4

%% Set up environment

addpath(genpath('bin/'))
colors = distinguishable_colors(101); colors(4,:) = [];

%% Load data

full_table = readtable('../data/Table_S2_embryo_data_table.csv');

%% REMOVE?

%bins = readtable('mm10_1Mb_bins.txt');
%gc = readtable('mm10_1Mb_gc_content.txt');

%% fix GC column REMOVE?

%for i=1:size(gc,1); disp(i)
%    
%    sel = (full_table.chr == gc.Var5(i) & full_table.pos >= gc.Var2(i) & full_table.pos < gc.Var3(i));
%    full_table.gc(sel) = gc.Var8(i);
%    
%end

%full_table.gc(full_table.gc>1) = 0;

%% Data selection

% select only inliers

sel_table = full_table(full_table.inlier == 1,:);

% create list of all embryos 

embryos = unique(sel_table.embryo_id)';
num_cells = size(unique(sel_table{:,[1 2]},'rows'),1);

%% Calculate separation score per read, chromosome, and cell

% set number of neighbors for analysis

k = 100;

% holds scores

separation_score = zeros(size(sel_table,1),1);
chr_separation_score = zeros(num_cells,21);
cell_separation_score = zeros(num_cells,1);

% holds cell stage

cell_stage = strings(num_cells,1);

%rabl_score = zeros(size(sel_table,1),1);

% loop through all embryos

cell_count = 1;
for embryo=[embryos]
    
    % find and loop through all cells in embryo
    
    embryo_table = sel_table(sel_table.embryo_id == embryo,:);
    cells = unique(embryo_table.cell_id)';
    
    for cell=[cells]
        
        % select data for cell
        
        cell_table = embryo_table(embryo_table.cell_id == cell,:);
        
        % find k nearest spatial neighbors
        
        knn_idx = knnsearch([cell_table.x_um_abs cell_table.y_um_abs cell_table.z_um_abs],[cell_table.x_um_abs cell_table.y_um_abs cell_table.z_um_abs],'K',k);
        
        % find neighbors on the same homolog for exclusion
        
        same_homolog = (cell_table.chr(knn_idx) == cell_table.chr) & (cell_table.cluster_hap_imputed(knn_idx) == cell_table.cluster_hap_imputed);
        
        % count maternal and paternal SNPs in neighbors
        
        mat_haps = sum(cell_table.cluster_hap_imputed(knn_idx) == 0 & (~same_homolog),2);
        pat_haps = sum(cell_table.cluster_hap_imputed(knn_idx) == 1 & (~same_homolog),2);
        
        % calculate separation scores
        
        separation_score(sel_table.embryo_id == embryo & sel_table.cell_id == cell) = max(mat_haps,pat_haps)./(mat_haps+pat_haps);
        cell_separation_score(cell_count) = nanmean(max(mat_haps,pat_haps)./(mat_haps+pat_haps));
        
        % remove separation scores for cells with low coverage
        
        if size(cell_table,1) < 500
            cell_separation_score(cell_count) = NaN;
        end
        
        % calculate separation score per chromosome
        
        for chr=1:21
            in_chr = (cell_table.chr == chr);
            chr_separation_score(cell_count,chr) = nanmean(max(mat_haps(in_chr),pat_haps(in_chr))./(mat_haps(in_chr)+pat_haps(in_chr)));
        end
            
        % set cell stage and increment cell count
        
        cell_stage(cell_count) = string(cell_table.stage(1));
        cell_count = cell_count + 1;
        
    end
end

%% Assign index to each stage

stage_idx = size(cell_stage,1);
stage_idx(string(cell_stage) == 'zygote') = 1;
stage_idx(string(cell_stage) == '2cell') = 2;
stage_idx(string(cell_stage) == '4cell') = 3;

%% Plot separation score by stage

figure; notBoxPlot(cell_separation_score,stage_idx);
xticklabels({'zygote','2-cell','4-cell'})
ylabel('Mean separation score per cell')
ylim([.7 1]);

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 4 4];

%% Perform K-S test for signficance values

[h,p1] = kstest2(cell_separation_score(stage_idx==1),cell_separation_score(stage_idx==2))
[h,p2] = kstest2(cell_separation_score(stage_idx==2),cell_separation_score(stage_idx==3))

%% Plot separation score by stage and chromosome

stages = ["zygote","2-cell","4-cell"];

figure; 

for i=1:3
    subplot(3,1,i)
    notBoxPlot(chr_separation_score(stage_idx==i,1:20));
    xticklabels(cat(1,string(1:19)',"X"))
    ylabel({'Mean separation score','per chromosome'})
    ylim([.7 1]); 
    text(10.5,1.05,char(stages(i)),'HorizontalAlignment','center','FontSize',14)
end

xlabel('Chromosome','FontSize',12)

%% Calculate Rabl score per read

cell_count = 1;
for embryo=[embryos]
    
    % find and loop through all cells in embryo
    
    embryo_table = sel_table(sel_table.embryo_id == embryo,:);
    cells = unique(embryo_table.cell_id)';
    
    for cell=[cells]
        
        % select data for cell
        
        cell_table = embryo_table(embryo_table.cell_id == cell,:);
        
        % find k nearest spatial neighbors
        
        knn_idx = knnsearch([cell_table.x_um_abs cell_table.y_um_abs cell_table.z_um_abs],[cell_table.x_um_abs cell_table.y_um_abs cell_table.z_um_abs],'K',k);
        
        % find neighbors on the same homolog for exclusion
        
        same_homolog = (cell_table.chr(knn_idx) == cell_table.chr) & (cell_table.cluster_hap_imputed(knn_idx) == cell_table.cluster_hap_imputed);
        
        % calculate Rabl score
        
        neighborhood_chr_pos = cell_table.rel_chr_pos(knn_idx); 
        neighborhood_chr_pos(neighborhood_chr_pos == 0 | same_homolog) = NaN;
        rabl_score(sel_table.embryo_id == embryo & sel_table.cell_id == cell) = nanmean(neighborhood_chr_pos,2);
            
        % set cell stage and increment cell count
        
        cell_stage(cell_count) = string(cell_table.stage(1));
        cell_count = cell_count + 1;
        
    end
end

%% Bin Rabl scores by chromosome position and stage

% set number of bins and window

num_bins = 100;
win = 0.5;

median_pos = zeros(num_bins,3);
for i=1:num_bins; disp(i)
    median_pos(i,1) = nanmedian(rabl_score(string(sel_table.stage)=='zygote' & sel_table.rel_chr_pos > (i-win)./100 & sel_table.rel_chr_pos < (i+win)./100));
    median_pos(i,2) = nanmedian(rabl_score(string(sel_table.stage)=='2cell' & sel_table.rel_chr_pos > (i-win)./100 & sel_table.rel_chr_pos < (i+win)./100));
    median_pos(i,3) = nanmedian(rabl_score(string(sel_table.stage)=='4cell' & sel_table.rel_chr_pos > (i-win)./100 & sel_table.rel_chr_pos < (i+win)./100));
end

%% Plot Rabl scores by stage

figure;
scatter(1:100,median_pos(:,1),'filled'); hold on;
scatter(1:100,median_pos(:,2),'filled'); hold on;
scatter(1:100,median_pos(:,3),'filled'); hold on;

ylim([.25 .75])
xticks([0:25:100])
yticks([0.25:0.25:0.75])

xticklabels({[0:0.25:1]})

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 4 4];

legend('zygote','2-cell','4-cell','Location','southeast')
xlabel('Centromere-telomere position')
ylabel({'Mean centromere-telomere';'position of spatial neighbors'})

%% Get correlations for each stage

corr(sel_table.rel_chr_pos(string(sel_table.stage)=='zygote'),rabl_score(string(sel_table.stage)=='zygote')')
corr(sel_table.rel_chr_pos(string(sel_table.stage)=='2cell'),rabl_score(string(sel_table.stage)=='2cell')')
corr(sel_table.rel_chr_pos(string(sel_table.stage)=='4cell'),rabl_score(string(sel_table.stage)=='4cell')')

%%


