% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Spatial separation, Rabl, and GC analyses for Fig. 4

%% Set up environment

addpath(genpath('bin/'))
colors = distinguishable_colors(101); colors(4,:) = [];

%% Load data

full_table = readtable('../data/Table_S2_embryo_data_table.csv');

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
        
        same_homolog = (cell_table.chr(knn_idx) == cell_table.chr) & (cell_table.cluster(knn_idx) == cell_table.cluster);
        
        % count maternal and paternal SNPs in neighbors
        
        mat_haps = sum(cell_table.hap_assignment(knn_idx) == 0 & (~same_homolog),2);
        pat_haps = sum(cell_table.hap_assignment(knn_idx) == 1 & (~same_homolog),2);
        
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

%% Get stats for main text

nanstd(cell_separation_score(stage_idx == 1))
nanstd(cell_separation_score(stage_idx == 2))
nanstd(cell_separation_score(stage_idx == 3))

[h,p1] = kstest2(cell_separation_score(stage_idx==1),cell_separation_score(stage_idx==2))
[h,p2] = kstest2(cell_separation_score(stage_idx==2),cell_separation_score(stage_idx==3))

%% Plot separation score by stage (Fig. 4C)

figure; notBoxPlot(cell_separation_score,stage_idx);
xticklabels({'zygote','2-cell','4-cell'})
ylabel('Mean separation score per cell')
ylim([.75 1]);

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 4 4];

%% Plot separation score by stage and chromosome (Fig. S21)

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
        
        same_homolog = (cell_table.chr(knn_idx) == cell_table.chr) & (cell_table.cluster(knn_idx) == cell_table.cluster);
        
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

% loop through bins and get median for each stage

median_pos = zeros(num_bins,3);
for i=1:num_bins; disp(i)
    median_pos(i,1) = nanmedian(rabl_score(string(sel_table.stage)=='zygote' & sel_table.rel_chr_pos > (i-win)./100 & sel_table.rel_chr_pos < (i+win)./100));
    median_pos(i,2) = nanmedian(rabl_score(string(sel_table.stage)=='2cell' & sel_table.rel_chr_pos > (i-win)./100 & sel_table.rel_chr_pos < (i+win)./100));
    median_pos(i,3) = nanmedian(rabl_score(string(sel_table.stage)=='4cell' & sel_table.rel_chr_pos > (i-win)./100 & sel_table.rel_chr_pos < (i+win)./100));
end

%% Plot Rabl scores by stage (Fig. 4E)

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

%% Create haplotype genomic bins

bins = readtable('ref/mm10_1Mb_bins.txt');
gc = readtable('ref/mm10_1Mb_gc.txt');

no_y = 1:min(find(string(bins.chr) == 'chrY'))-1;
hap_bins = [bins(no_y,:); bins(no_y,:)];

hap_bins.hap = zeros(size(hap_bins,1),1);
hap_bins.hap(1:size(no_y,2)) = 0;
hap_bins.hap(size(no_y,2)+1:end) = 1;

hap_bins.bin_ind = (1:size(hap_bins,1))';

hap_bins.gc = [gc.Var8(no_y); gc.Var8(no_y)];
hap_bins.gc(hap_bins.gc == 0) = NaN;

%% Assign zygote reads to a haplotyped bin

% get list of zygotes

zygote_table = sel_table(string(sel_table.stage) == 'zygote',:);
zygotes = [unique(zygote_table.embryo_id)]';
num_zygotes = size(unique([zygote_table.embryo_id zygote_table.cell_id],'rows'),1);

% assign reads

zygote_table.hap_bin = zeros(size(zygote_table,1),1);

for i=1:size(zygote_table,1);
    bin_ind = hap_bins.bin_ind(hap_bins.chr_ind == zygote_table.chr(i) & hap_bins.bin_start < zygote_table.pos(i) ...
        & hap_bins.bin_end > zygote_table.pos(i) & hap_bins.hap == zygote_table.cluster_hap_imputed(i));
    if size(bin_ind,1) == 1
        zygote_table.hap_bin(i) = bin_ind;
    end
end

%% Create zygote distance matrices for lamina and NPBs

zygote_lam_dist_mat = NaN([size(hap_bins,1) num_zygotes]);
zygote_npb_dist_mat = NaN([size(hap_bins,1) num_zygotes]);

% loop through all zygotes

cell_count = 1;
for embryo=[zygotes]; disp(embryo)
    
    embryo_table = zygote_table(zygote_table.embryo_id == embryo,:);
    cells = unique(embryo_table.cell_id)';
    
    for cell=[cells]
        
        % select all haplotyped chromosomes
        
        cell_table = embryo_table(embryo_table.cell_id == cell & (embryo_table.cluster_hap_imputed == 0 | embryo_table.cluster_hap_imputed == 1) & embryo_table.chr <= 20,:);
        max_bin = max(cell_table.hap_bin(:));

        % fill in distance matrices
        
        zygote_lam_dist_mat(1:max_bin,cell_count) = accumarray(cell_table.hap_bin(:),cell_table.dist_to_lamin)./accumarray(cell_table.hap_bin(:),1);
        zygote_npb_dist_mat(1:max_bin,cell_count) = accumarray(cell_table.hap_bin(:),cell_table.dist_to_npb)./accumarray(cell_table.hap_bin(:),1);

        cell_count = cell_count+1;

    end     
    
end

%% Assign 2-cell reads to a haplotyped bin

% get list of 2-cells

twocell_table = sel_table(string(sel_table.stage) == '2cell',:);
twocells = [unique(twocell_table.embryo_id)]';
num_twocells = size(unique([twocell_table.embryo_id twocell_table.cell_id],'rows'),1);

% assign reads

twocell_table.hap_bin = zeros(size(twocell_table,1),1);

for i=1:size(twocell_table,1); disp(i)
    bin_ind = hap_bins.bin_ind(hap_bins.chr_ind == twocell_table.chr(i) & hap_bins.bin_start < twocell_table.pos(i) ...
        & hap_bins.bin_end > twocell_table.pos(i) & hap_bins.hap == twocell_table.cluster_hap_imputed(i));
    if size(bin_ind,1) == 1
        twocell_table.hap_bin(i) = bin_ind;
    end
end

%% Create 2-cell distance matrices for lamina and NPBs

twocell_lam_dist_mat = NaN([size(hap_bins,1) num_twocells]);
twocell_npb_dist_mat = NaN([size(hap_bins,1) num_twocells]);

% loop through all 2-cells

cell_count = 1;
for embryo=[twocells]; disp(embryo)
    
    embryo_table = twocell_table(twocell_table.embryo_id == embryo,:);
    cells = unique(embryo_table.cell_id)';
    
    for cell=[cells]
        
        % select all haplotyped chromosomes
        
        cell_table = embryo_table(embryo_table.cell_id == cell & (embryo_table.cluster_hap_imputed == 0 | embryo_table.cluster_hap_imputed == 1) & embryo_table.chr <= 20,:);
        max_bin = max(cell_table.hap_bin(:));

        % fill in distance matrices
        
        twocell_lam_dist_mat(1:max_bin,cell_count) = accumarray(cell_table.hap_bin(:),cell_table.dist_to_lamin)./accumarray(cell_table.hap_bin(:),1);
        twocell_npb_dist_mat(1:max_bin,cell_count) = accumarray(cell_table.hap_bin(:),cell_table.dist_to_npb)./accumarray(cell_table.hap_bin(:),1);

        cell_count = cell_count+1;

    end     
    
end

%% Plot correlation of GC and lamina for maternal/paternal chr12 (Fig. 4G)

chr = 12;
chr_size = size(find(bins.chr_ind == chr),1)

% calculate correlations

mat_corr = corr(nanmean(zygote_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 0,:),2),hap_bins.gc(hap_bins.chr_ind == chr & hap_bins.hap == 0),'type','Spearman','rows','complete');
pat_corr = corr(nanmean(zygote_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 1,:),2),hap_bins.gc(hap_bins.chr_ind == chr & hap_bins.hap == 1),'type','Spearman','rows','complete');

% plot maternal

figure;
plot(1:chr_size,nanmean(zygote_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 0,:),2),'Color',[147 39 143]./255,'LineWidth',2); hold on;
xlim([0 121]); ylim([0 4]); yticks(0:1:4); ylabel({'Mean distance','to lamina (\mum)'})

yyaxis right
plot(1:chr_size,gc.Var8(gc.Var5 == chr),'Color',[217 83 25]./255,'LineWidth',2); hold on;
ylim([.35 .55]); ylabel('GC content')

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 8 2];

% plot paternal

figure;
plot(1:chr_size,nanmean(zygote_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 1,:),2),'Color',[0 146 69]./255,'LineWidth',2); hold on;
xlim([0 121]); ylim([0 4]); yticks(0:1:4); ylabel({'Mean distance','to lamina (\mum)'})

yyaxis right
plot(1:chr_size,gc.Var8(gc.Var5 == chr),'Color',[217 83 25]./255,'LineWidth',2); hold on;
ylim([.35 .55]); ylabel('GC content')

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 8 2];

%% Calculate correlation between GC content and distance to nuclear lamina per chromosome

chr_corrs = zeros(20,4);

for chr=1:20
    
    x = gc.Var8(find(hap_bins.chr_ind == chr & hap_bins.hap == 0)); y = nanmean(zygote_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 0,:),2);
    chr_corrs(chr,1) = corr(x(~isnan(x) & ~isnan(y)),y(~isnan(x) & ~isnan(y)),'type','Spearman');

    x = gc.Var8(find(hap_bins.chr_ind == chr & hap_bins.hap == 0)); y = nanmean(zygote_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 1,:),2);
    chr_corrs(chr,2) = corr(x(~isnan(x) & ~isnan(y)),y(~isnan(x) & ~isnan(y)),'type','Spearman');
    
    x = gc.Var8(find(hap_bins.chr_ind == chr & hap_bins.hap == 0)); y = nanmean(twocell_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 0,:),2);
    chr_corrs(chr,3) = corr(x(~isnan(x) & ~isnan(y)),y(~isnan(x) & ~isnan(y)),'type','Spearman');

    x = gc.Var8(find(hap_bins.chr_ind == chr & hap_bins.hap == 0)); y = nanmean(twocell_lam_dist_mat(hap_bins.chr_ind == chr & hap_bins.hap == 1,:),2);
    chr_corrs(chr,4) = corr(x(~isnan(x) & ~isnan(y)),y(~isnan(x) & ~isnan(y)),'type','Spearman');
    
end

%% Plot GC-distance to lamina correlation by stage and haplotype (Fig. 4H)

figure; 
h= notBoxPlot(chr_corrs);
xticklabels({'mat zygote','pat zygote','mat 2-cell','pat 2-cell'})
d=[h.data];
set(d([1 3]),'markerfacecolor',[147 39 143]./255,'color','black')
set(d([2 4]),'markerfacecolor',[0 146 69]./255,'color','black')

ylim([0 1])
[h_z p_z] = kstest2(chr_corrs(:,1),chr_corrs(:,2));
[h_t p_t] = kstest2(chr_corrs(:,3),chr_corrs(:,4));

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 4 4];

%% Plot genomic distance to lamaina/NPBs for zygotes and 2-cells (Fig. S22)

figure; 

% zygote distance to lamina
scatter(no_y,nanmean(zygote_lam_dist_mat(hap_bins.hap==0,:),2),[],[147 39 143]./255,'.'); hold on;
scatter(no_y,nanmean(zygote_lam_dist_mat(hap_bins.hap==1,:),2),[],[0 146 69]./255,'.'); hold on;

% zygote distance to NPBs
%scatter(no_y,nanmean(zygote_npb_dist_mat(hap_bins.hap==0,:),2),[],[147 39 143]./255,'.'); hold on;
%scatter(no_y,nanmean(zygote_npb_dist_mat(hap_bins.hap==1,:),2),[],[0 146 69]./255,'.'); hold on;

% 2-cell distance to lamina
%scatter(no_y,nanmean(twocell_lam_dist_mat(hap_bins.hap==0,:),2),[],[147 39 143]./255,'.'); hold on;
%scatter(no_y,nanmean(twocell_lam_dist_mat(hap_bins.hap==1,:),2),[],[0 146 69]./255,'.'); hold on;

% 2-cell distance to NPBs
%scatter(no_y,nanmean(twocell_npb_dist_mat(hap_bins.hap==0,:),2),[],[147 39 143]./255,'.'); hold on;
%scatter(no_y,nanmean(twocell_npb_dist_mat(hap_bins.hap==1,:),2),[],[0 146 69]./255,'.'); hold on;

xlim([0 size(zygote_lam_dist_mat(hap_bins.hap==0,:),1)])
ylim([0 10])

text_y = ylims(2) + ylims(2)/15;
tick_y = [ylims(2) - ylims(2)/50, ylims(2) + ylims(2)/50];

bin_size = 1000000;
num_chrs = max(bins.chr_ind(no_y));

chr_lens = zeros(num_chrs,1);
for i=1:num_chrs
   chr_lens(i) = max(bins.bin_end(bins.chr_ind == i));
end

count = 1;
text_col = 'black';
for chr=1:num_chrs

    count = count+ceil((chr_lens(chr)/bin_size));
    
    plot(repmat([count],2,1),tick_y,'Color','black','Marker','none','HandleVisibility','off'); hold on;
    if chr == 1
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'chr1','Color',text_col,'HorizontalAlignment','center','FontSize',12); hold on;
    elseif chr == num_chrs
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'X','Color',text_col,'HorizontalAlignment','center','FontSize',12); hold on;
    elseif chr >= 0 & chr <= 8
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,sprintf('%d',chr),'Color',text_col,'HorizontalAlignment','center','FontSize',12); hold on;
    elseif chr == 9
        t = text(count-(ceil(chr_lens(chr)/bin_size)./2),text_y,'...','Color',text_col,'HorizontalAlignment','center','FontSize',12); hold on;
    end
    
end

xticks([])
a = get(gca,'YTickLabel');  set(gca,'fontsize',12)
xlabel('Genomic position','FontSize',12); 

ylabel({'Distance to nuclear', 'lamina (\mum)'},'FontSize',12)
%ylabel({'Distance to nucleolus','precursor bodies (\mum)'},'FontSize',12)

fig = gcf;
fig.Units = 'inches';
fig.Position = [0 0 6 2];  
