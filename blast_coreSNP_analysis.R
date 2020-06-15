setwd("/Volumes/data/MPV/THEJ/Projekter/Staphylococcus/Staphylococcus_epidermidis_point_mutation_resistance/")
library(reshape2)
library(pls)
library(ggplot2)

SNP_table = read.table("nasp_assemblies_1162/matrices/missingdata_SNP.tsv",sep="\t",header=TRUE,check.names = F)
pa_table = read.table("nasp_assemblies_1162/matrices/missingdata_presence.tsv",sep="\t",header=TRUE,check.names = F)

SNP_table$contig = unlist(lapply(as.vector(SNP_table$LocusID),function(x) strsplit(x,'::')[[1]][1]))
SNP_table$pos = as.numeric(unlist(lapply(as.vector(SNP_table$LocusID),function(x) strsplit(x,'::')[[1]][2])))
pa_table$contig = unlist(lapply(as.vector(pa_table$LocusID),function(x) strsplit(x,'::')[[1]][1]))
pa_table$pos = as.numeric(unlist(lapply(as.vector(pa_table$LocusID),function(x) strsplit(x,'::')[[1]][2])))


master_SNP[1:10,1:10]
gene_table[1:10,]

gene_table = read.table("ATCC12228_cds_presence_tab.txt",sep = "\t",header=TRUE)

snp_pca = prcomp()


gene_table$var = 1
gene_presence_df = dcast(gene_table,formula = gene~ID,value.var = 'var')
rownames(gene_presence_df) = gene_presence_df$gene
gene_presence_mat = as.matrix(gene_presence_df[,-1])
gene_presence_mat[gene_presence_mat>0] = 1
gene_presence_df = as.data.frame(gene_presence_mat)
gene_presence_df[1:10,1:10]


pca = prcomp(t(gene_presence_df))
plot_df = data.frame(pca$x)

labels = c(paste0("PC1 [",sprintf("%.1f",explvar(pca)[1]),"%]"),paste0("PC2 [",sprintf("%.1f",explvar(pca)[2]),"%]"))

p = ggplot(as.data.frame(plot_df),aes(x=PC1,y=PC2,color = type)) + labs(color = "gene presence") + geom_point(size=1.5, alpha=1)+ stat_ellipse(level=0.75) + xlab(labels[1]) + ylab(labels[2]) + theme_bw()
p
p = ggplot(as.data.frame(plot_df),aes(x=PC3,y=PC4,color = type)) + labs(color = "gene presence") + geom_point(size=1.5, alpha=1)+ stat_ellipse(level=0.75) + xlab(labels[1]) + ylab(labels[2]) + theme_bw()
p

grep('NZ_CP022247.1_cds_WP_002447641.1_1857',rownames(gene_presence_df))

grep('NZ_CP022247.1_cds_1804__complement_1870445',rownames(gene_presence_df))




load_SNP_table = function


get_SNP_dist <- function(filename,SNP_table,pa_table) {
  namesplit = strsplit(filename,'__')[[1]]
  posstring = namesplit[2]
  contig = paste0(strsplit(filename,'_')[[1]][1:2],collapse='_')
  if (substr(posstring,1,4) != "join" & substr(posstring,1,15) != "complement_join") {
    rgx = regexec("[complement_]*[A-Z]*([0-9]*)\\.\\.[A-Z]*([0-9]*)_*\\.txt",posstring)
    rgx_matches = regmatches(posstring,rgx)[[1]]
    startpos = as.numeric(rgx_matches[[2]])
    endpos = as.numeric(rgx_matches[[3]])
    if (startpos > endpos) {
      tmp = startpos
      startpos = endpos
      endpos = tmp
    }
    startpos = startpos-10000
    endpos = endpos+10000
    if (startpos<1) {
      startpos = 1
    }
    include_idx = which(SNP_table$contig==contig & SNP_table$pos>startpos & SNP_table$pos<endpos)
    sub_SNP_table = SNP_table[include_idx,-c(1,2,(ncol(SNP_table)-1),ncol(SNP_table))]
    d = dist(t(sub_SNP_table),method = "binary")
    return_obj = list('check'=1,'dist_obj'=d,'snp_count'=length(include_idx),"sub_snp_table"=sub_SNP_table)
  }
  else {
    return_obj= list('check'=0)
  }
  return(return_obj)
}

get_pairwise_distance <- function(dist_mat,ID1, ID2 = NA) {
  if (is.na(ID2)) {
    sub_dist_mat = dist_mat[which(rownames(dist_mat) %in% ID1),which(colnames(dist_mat) %in% ID1)]
    dist_vec = sub_dist_mat[upper.tri(sub_dist_mat)]
  } else {
    dist_vec = as.vector(dist_mat[which(rownames(dist_mat) %in% ID1),which(colnames(dist_mat) %in% ID2)])
  }
  return(dist_vec)
}

load_nonsym_SNP_table <- function(filename) {
  tbl = read.table(filename,header=TRUE,sep = "\t")
  return(tbl)
}

get_nonsym_SNP_pa_position <- function(IDs,tbl,position) {
  return_vec = rep("No SNP",length(IDs))
  with_SNP_IDs = as.vector(tbl$ID)[which(tbl$var_pos==position)]
  return_vec[which(IDs %in% with_SNP_IDs)] = "SNP"
  return_factor = factor(return_vec,levels=c("SNP","No SNP"))
  return(return_factor)
}

get_nonsym_SNP_pa_variant <- function(IDs,tbl,variant) {
  return_vec = rep("No SNP",length(IDs))
  with_SNP_IDs = as.vector(tbl$ID)[which(tbl$var==variant)]
  return_vec[which(IDs %in% with_SNP_IDs)] = "SNP"
  return_factor = factor(return_vec,levels=c("SNP","No SNP"))
  return(return_factor)
}

SNP_boxplot <- function(dist_mat,SNP_factor) {
  IDs = rownames(dist_mat)
  pos_IDs = IDs[which(SNP_factor=="SNP")]
  neg_IDs = IDs[which(SNP_factor=="No SNP")]
  d_pos = get_pairwise_distance(dist_mat,pos_IDs)
  d_neg = get_pairwise_distance(dist_mat,neg_IDs)
  d_pn = get_pairwise_distance(dist_mat,pos_IDs,neg_IDs)
  plot_dist_table = data.frame("Distance"=c(d_pos,d_neg,d_pn),"Type"=c(rep("Within SNP isolates",length(d_pos)),rep("Within no SNP isolates",length(d_neg)),rep("Between SNP/no SNP isolates",length(d_pn))))
  p <- ggplot(plot_dist_table,aes(x=Type,y=Distance,color=Type)) + geom_boxplot() + scale_color_manual(values = RColorBrewer::brewer.pal(3,"Set1")) + geom_jitter(mapping = NULL, data = NULL, stat = "identity",position = "jitter")
  return(list('plot'=p,'plot_table'=plot_dist_table))
}

SNP_histogram <- function(dist_mat,SNP_factor) {
  IDs = rownames(dist_mat)
  pos_IDs = IDs[which(SNP_factor=="SNP")]
  neg_IDs = IDs[which(SNP_factor=="No SNP")]
  d_pos = get_pairwise_distance(dist_mat,pos_IDs)
  d_neg = get_pairwise_distance(dist_mat,neg_IDs)
  d_pn = get_pairwise_distance(dist_mat,pos_IDs,neg_IDs)
  plot_dist_table = data.frame("Distance"=c(d_pos,d_neg,d_pn),"Type"=c(rep("Within SNP isolates",length(d_pos)),rep("Within no SNP isolates",length(d_neg)),rep("Between SNP/no SNP isolates",length(d_pn))))
  p <- ggplot(plot_dist_table,aes(x=Distance,fill=Type)) + geom_histogram(position="dodge") + scale_fill_manual(values = RColorBrewer::brewer.pal(3,"Set1"))
  return(list('plot'=p,'plot_table'=plot_dist_table))
}

SNP_index <- function(dist_mat,SNP_factor) {
  IDs = rownames(dist_mat)
  pos_IDs = IDs[which(SNP_factor=="SNP")]
  neg_IDs = IDs[which(SNP_factor=="No SNP")]
  d_pos = get_pairwise_distance(dist_mat,pos_IDs)
  d_neg = get_pairwise_distance(dist_mat,neg_IDs)
  d_pn = get_pairwise_distance(dist_mat,pos_IDs,neg_IDs)
  rank_order = order(c(d_pos,d_neg,d_pn))
  r_pos = rank_order[1:length(d_pos)]
  r_neg = rank_order[(length(d_pos)+1):(length(d_pos)+length(d_neg))]
  r_pn = rank_order[(length(d_pos)+length(d_neg)+1):length(rank_order)]
  return(c(mean(d_pos),mean(d_neg),mean(c(d_pos,d_neg)),mean(d_pn),mean(r_pos),mean(r_neg),mean(c(r_pos,r_neg)),mean(r_pn)))
  #return(data.frame('Distances'=c(d_pos,d_neg,d_pn), 'ranks' = c(r_pos,r_neg,r_pn),"Type"=c(rep("Within SNP isolates",length(d_pos)),rep("Within no SNP isolates",length(d_neg)),rep("Between SNP/no SNP isolates",length(d_pn)))))
}

add_to_df_and_print_itol <- function(d,tbl,pos,refAA,AA,gene_name) {
  mut = paste0(refAA,pos,AA)
  add_vec = rep(refAA,nrow(d))
  snp_IDs = as.vector(tbl$ID)[which(tbl$var==mut)]
  add_vec[which(d$Sequence.Id.used.in.this.study %in% snp_IDs)] = AA
  add_vec = factor(add_vec)
  d$add_vec = add_vec
  gene_mut = paste0(gene_name,'_',mut)
  colnames(d)[ncol(d)] = gene_mut
  itol_lines = setup_color_lines(d$Sequence.Id.used.in.this.study,add_vec,RColorBrewer::brewer.pal(2,"Set1"),gene_mut)
  print_template(template,itol_lines,paste0("international_",gene_mut,"_COLORSTRIP.txt"),legend_title = gene_mut)
  return(d)
}

#### Set up SNP PA matrix ####
files = list.files("ATCC12228_cds_SNP_summaries")
pa_mat = matrix(nrow=0,ncol=7)
for (i in 1:length(files)) {
  filename = files[i]
  print(i)
  if (filename %in% rownames(gene_presence_df)) {
    presence_sum = sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
    print(presence_sum)
    if (presence_sum>1140) {
      tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
      if (nrow(tbl)>0) {
        tbl$file = filename
        tbl = as.matrix(tbl)
        pa_mat = rbind(pa_mat,tbl)
      }
      
    }
  }
  
  
}

pa_tbl = as.data.frame(pa_mat)
pa_tbl$value = 1
pa_tbl$gene_pos = paste0(as.vector(pa_tbl$file),'::',pa_tbl$var_pos)
SNP_pa_df = dcast(pa_tbl,formula = gene_pos~ID,value.var = 'value')
rownames(SNP_pa_df) = SNP_pa_df$gene_pos
snp_pa_mat = as.matrix(SNP_pa_df[,-1])
snp_pa_mat[snp_pa_mat>0] = 1
snp_pa_df = as.data.frame(snp_pa_mat)

write.table(snp_pa_df,"nonsym_SNP_presence_matrix.txt",sep = "\t")

##################

files = list.files("ATCC12228_cds_SNP_summaries")
stats_mat = matrix(nrow=0,ncol=13)
for (i in 1:length(files)) {
  filename = files[i]
  tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
  snp_pos_prevalence = sort(table(tbl$var_pos))
  pos_to_check = names(snp_pos_prevalence)[which(snp_pos_prevalence > 50 & snp_pos_prevalence < 1111)]
  if (length(pos_to_check) > 0) {
    snp_dist_obj = get_SNP_dist(filename,SNP_table = SNP_table,pa_table = pa_table)
    if (snp_dist_obj$check == 1) {
      dist_obj = snp_dist_obj$dist_obj
      snp_count = snp_dist_obj$snp_count
      dist_mat = as.matrix(dist_obj)
      print(paste0("parsing gene ",i," ",filename," ",length(pos_to_check)," SNP positions")) %>% 
      for (pos in pos_to_check) {
        SNP_pa_factor = get_nonsym_SNP_pa_position(rownames(dist_mat),tbl,pos)
        pos_count = length(which(SNP_pa_factor=="SNP"))
        neg_count = length(which(SNP_pa_factor=="No SNP"))
        snp_stats = SNP_index(dist_mat,SNP_pa_factor)
        add_vec = c(filename,pos,snp_count,pos_count,neg_count,snp_stats)
        print(length(add_vec))
        stats_mat = rbind(stats_mat,add_vec)
      }
    }
  } else {
    print(paste0("Skipping ",filename," no positions with SNPS within threshold prevalence"))
  }
  
}
stats_df = as.data.frame(stats_mat)
colnames(stats_df) = c("file","pos","SNPs_in_range","Isolates_with_SNP","Isolates_without_SNP","d_pos","d_neg","d_within","d_pn","r_pos","r_neg","r_within","r_pn")
stats_df$r_neg = as.numeric(as.vector(stats_df$r_neg))
stats_df$r_pos = as.numeric(as.vector(stats_df$r_pos))
stats_df$r_pn = as.numeric(as.vector(stats_df$r_pn))
stats_df$d_neg = as.numeric(as.vector(stats_df$d_neg))
stats_df$d_pos = as.numeric(as.vector(stats_df$d_pos))
stats_df$d_pn = as.numeric(as.vector(stats_df$d_pn))
stats_df$r_within = as.numeric(as.vector(stats_df$r_within))
stats_df$d_within = as.numeric(as.vector(stats_df$d_within))


#### setup SNP presence data frame ####
files = list.files("ATCC12228_cds_SNP_summaries")
nonsym_snp_mat = matrix(nrow=0,ncol=length(IDs))
colnames(nonsym_snp_mat) = IDs
gene_pos_vec = c()
for (i in 1:length(files)) {
  filename = files[i]
  tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
  snp_pos_prevalence = sort(table(tbl$var_pos))
  pos_to_check = names(snp_pos_prevalence)[which(snp_pos_prevalence > 50 & snp_pos_prevalence < 1111)]
  if (length(pos_to_check) > 0) {
    for (pos in pos_to_check) {
      gene_pos = paste0(filename,'__',pos)
      gene_pos_vec = c(gene_pos_vec,gene_pos)
      IDs_with_SNP = as.vector(tbl$ID)[which(tbl$var_pos==pos)]
      add_vec = rep(0,ncol(nonsym_snp_mat))
      add_vec[which(colnames(nonsym_snp_mat) %in% IDs_with_SNP)] = 1
      nonsym_snp_mat = rbind(nonsym_snp_mat,add_vec)
    }
  } else {
    print(paste0("Skipping ",filename," no positions with SNPS within threshold prevalence"))
  }
}


p = ggplot(stats_df,aes(x=d_pos,d_neg,color=file)) + geom_point() + theme(legend.position = "none")
p

stats_df[which(stats_df$r_pos/stats_df$r_neg>1),]
stats_df[which(stats_df$d_pos/mean(c(stats_df$d_neg,stats_df$d_pos,stats_df$d_pn))>1),]
stats_df[which(stats_df$d_pos/stats_df$d_pn>0.9),]


stats_df = read.table("nonsym_SNPS_distances.txt",sep = "\t")
stats_df[1:10,]


stats_df$d_score = stats_df$d_pn/stats_df$d_within
stats_df$r_score = stats_df$r_pn/stats_df$r_within


head(stats_df[order(stats_df$d_score,decreasing = T),])

filename = "NZ_CP022247.1_cds_WP_088922884.1_1902__1972649..1974670.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==643),]
d_new = add_to_df_and_print_itol(d_new,tbl,643,"I","V",gene_name = "copper-translocating")


mean_d_score = unlist(lapply(levels(stats_df$file),function(x) mean(stats_df$d_score[which(stats_df$file==x)])))
names(mean_d_score) = levels(stats_df$file)
stats_df$mean_gene_d_score = NA
for (ID in names(mean_d_score)) {
  stats_df$mean_gene_d_score[which(stats_df$file==ID)] = mean_d_score[which(names(mean_d_score)==ID)]
}

stats_df$gene_adjusted_d_score = stats_df$d_score/stats_df$mean_gene_d_score


mean_rank_score = unlist(lapply(levels(stats_df$file),function(x) mean(stats_df$r_score[which(stats_df$file==x)])))
names(mean_rank_score) = levels(stats_df$file)
stats_df$mean_gene_r_score = NA
for (ID in names(mean_rank_score)) {
  stats_df$mean_gene_r_score[which(stats_df$file==ID)] = mean_rank_score[which(names(mean_rank_score)==ID)]
}

stats_df$gene_adjusted_r_score = stats_df$r_score/stats_df$mean_gene_r_score

gene_headers = readLines("reference_ATCC12228/GCF_002215535.1_ASM221553v1_cds_from_genomic_headers.txt")
header_split_1 = unlist(lapply(gene_headers, function(x) strsplit(x,'|',fixed = T)[[1]][2]))
header_split_2 = unlist(lapply(header_split_1, function(x) strsplit(x,' ',fixed = T)[[1]][1]))
names(header_split_2) = header_split_1
stats_df$gene = unlist(lapply(as.vector(stats_df$file), function(x) strsplit(x,'__',fixed=T)[[1]][1]))

idx_list = match(stats_df$gene,header_split_2)
stats_df$fasta_header = header_split_1[idx_list]

write.table(stats_df,file = "nonsym_SNP_statistics.txt",sep = "\t")

plot(stats_df$r_score,stats_df$gene_adjusted_r_score)
plot(stats_df$d_score,stats_df$gene_adjusted_d_score)

plot(stats_df$d_score,stats_df$r_score)


stats_df[which(stats_df$d_score>30),]
stats_df[which(stats_df$gene_adjusted_d_score>3),]

stats_df[which(stats_df$gene_adjusted_r_score>2),]

stats_df[which(stats_df$r_score>1.2),]

plot(stats_df$r_score,stats_df$gene_adjusted_r_score)
stats_df[which(stats_df$r_score>1.2),]

sub_df = stats_df[which(stats_df$Isolates_with_SNP>300 & stats_df$Isolates_without_SNP>300),]

plot(sub_df$r_score,sub_df$gene_adjusted_r_score)

sub_df[which(sub_df$d_score>5),]
sub_df[which(sub_df$gene_adjusted_d_score>1.5),]

sub_df[which(sub_df$d_score>2.5),]
sub_df[which(sub_df$gene_adjusted_d_score>1.5),]

d_new = add_to_df_and_print_itol(d_new,tbl,122,"A","G",gene_name = "DUF805")

d_new=d

filename = "NZ_CP022247.1_cds_2365__2495724..2496956.txt"
pos = 125
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==pos),]
d_new = add_to_df_and_print_itol(d_new,tbl,pos,"T","X",gene_name = "copper_resistance_protein")


filename = "NZ_CP022247.1_cds_WP_023474120.1_2224__2339432..2340295.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==153),]
d_new = add_to_df_and_print_itol(d_new,tbl,153,"I","M",gene_name = "yitT")

filename = "NZ_CP022247.1_cds_WP_088922733.1_1074__complement_1127428..1129776_.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==141),]
d_new = add_to_df_and_print_itol(d_new,tbl,141,"T","K",gene_name = "Serine_protease")

filename = "NZ_CP022247.1_cds_WP_088922733.1_1074__complement_1127428..1129776_.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==141),]
d_new = add_to_df_and_print_itol(d_new,tbl,141,"T","K",gene_name = "Serine_protease")

filename = "NZ_CP022247.1_cds_WP_002467876.1_1248__1301268..1301633.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==98),]
d_new = add_to_df_and_print_itol(d_new,tbl,98,"Q","K",gene_name = "globin")

filename = "NZ_CP022247.1_cds_WP_088922858.1_2280__2397282..2398028.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==26),]
d_new = add_to_df_and_print_itol(d_new,tbl,26,"E","A",gene_name = "glycosyltransferase")
tbl[which(tbl$var_pos==61),]
d_new = add_to_df_and_print_itol(d_new,tbl,61,"A","V",gene_name = "glycosyltransferase")

filename = "NZ_CP022247.1_cds_WP_002446111.1_1198__complement_1260257..1261012_.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==94),]
d_new = add_to_df_and_print_itol(d_new,tbl,94,"R","S",gene_name = "esterase")

filename = "NZ_CP022247.1_cds_WP_088922733.1_1074__complement_1127428..1129776_.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
tbl[which(tbl$var_pos==141),]
d_new = add_to_df_and_print_itol(d_new,tbl,141,"T","K",gene_name = "mutS2")

stats_df[which(stats_df$d_pos>0.8),]


plot(stats_df$r_score,stats_df$Isolates_with_SNP)
plot(stats_df$r_score,stats_df$gene_adjusted_r_score)
stats_df[which(stats_df$d_score>0.45 & stats_df$mean_),]

stats_df[which(stats_df$r_score>1.2 & stats_df$gene_adjusted_r_score>1),]


filename = "NZ_CP022247.1_cds_WP_002446111.1_1198__complement_1260257..1261012_.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])

tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
dist_obj = get_SNP_dist(filename,SNP_table = SNP_table,pa_table = pa_table)
sub_snp_table = dist_obj$sub_snp_table
dist_mat = as.matrix(dist_obj$dist_obj)
pa = get_nonsym_SNP_pa_position(IDs,tbl,94)

p = SNP_histogram(dist_mat,pa)
p
p = SNP_boxplot(dist_mat,pa)
p


pa = get_nonsym_SNP_pa_position(IDs,tbl,139)

p = SNP_histogram(dist_mat,pa)
p
p = SNP_boxplot(dist_mat,pa)
p

plot(rowMeans(sub_snp_table[,which(colnames(sub_snp_table) %in% snp_IDs)]))
plot(rowMeans(sub_snp_table[,which(!colnames(sub_snp_table) %in% snp_IDs)]))


filename = "NZ_CP022247.1_cds_WP_002447641.1_1857__1922352..1925033.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])

tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
dist_obj = get_SNP_dist(filename,SNP_table = SNP_table,pa_table = pa_table)
sub_snp_table = dist_obj$sub_snp_table
dist_mat = as.matrix(dist_obj$dist_obj)
pa = get_nonsym_SNP_pa_position(rownames(dist_mat),tbl,894)

p = SNP_histogram(dist_mat,pa)
p$plot
p = SNP_boxplot(dist_mat,pa)
p$plot


filename = "NZ_CP022247.1_cds_WP_002447895.1_1681__complement_1749002..1749820_.txt"
sum(gene_presence_df[which(rownames(gene_presence_df)==filename),])

tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
dist_obj = get_SNP_dist(filename,SNP_table = SNP_table,pa_table = pa_table)
sub_snp_table = dist_obj$sub_snp_table
dist_mat = as.matrix(dist_obj$dist_obj)
pa = get_nonsym_SNP_pa_position(IDs,tbl,184)

p = SNP_histogram(dist_mat,pa)
p
p = SNP_boxplot(dist_mat,pa)
p



filename = "ATCC12228_cds_SNP_summaries/NZ_CP022247.1_cds_WP_002447641.1_1857__1922352..1925033.txt"
tbl = load_nonsym_SNP_table(filename)

tbl = load_nonsym_SNP_table("ATCC12228_cds_SNP_summaries/NZ_CP022247.1_cds_WP_002447641.1_1857__1922352..1925033.txt")
dist_obj = get_SNP_dist("NZ_CP022247.1_cds_WP_001829751.1_158__139583..140020.txt",SNP_table = SNP_table,pa_table = pa_table)
dist_mat = as.matrix(dist_obj$dist_obj)
pa = get_nonsym_SNP_pa_position(IDs,tbl,76)

p = SNP_histogram(dist_mat,pa)
p
p = SNP_boxplot(dist_mat,pa)
p

tbl = load_nonsym_SNP_table("ATCC12228_cds_SNP_summaries/NZ_CP022247.1_cds_WP_002457296.1_2__270..686.txt")

dist_obj = get_SNP_dist("NZ_CP022247.1_cds_WP_002457296.1_2__270..686.txt",SNP_table = SNP_table,pa_table = pa_table)
dist_mat = as.matrix(dist_obj)
pa = get_nonsym_SNP_pa(IDs,tbl,43)

p = SNP_boxplot(dist_mat,pa)
p


filename = "NZ_CP022247.1_cds_WP_002495668.1_309__complement_297232..297795_.txt"


filename = "NZ_CP022247.1_cds_WP_002457296.1_2__270..686.txt"

tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
dist_obj = get_SNP_dist(filename,SNP_table = SNP_table,pa_table = pa_table)
dist_mat = as.matrix(dist_obj$dist_obj)
pa = get_nonsym_SNP_pa_position(IDs,tbl,25)

p = SNP_histogram(dist_mat,pa)
p
p = SNP_boxplot(dist_mat,pa)
p



### gyrA
filename = "NZ_CP022247.1_cds_WP_002447641.1_1857__1922352..1925033.txt"
tbl = load_nonsym_SNP_table(paste0("ATCC12228_cds_SNP_summaries/",filename))
dist_obj = get_SNP_dist(filename,SNP_table = SNP_table,pa_table = pa_table)
dist_mat = as.matrix(dist_obj$dist_obj)
pa = get_nonsym_SNP_pa_position(IDs,tbl,855)

p = SNP_histogram(dist_mat,pa)
p$plot
p = SNP_boxplot(dist_mat,pa)
p$plot


pa = get_nonsym_SNP_pa(IDs,tbl,548)
p = SNP_histogram(dist_mat,pa)
p
p = SNP_boxplot(dist_mat,pa)


pa = get_nonsym_SNP_pa(IDs,tbl,84)
p = SNP_boxplot(dist_mat,pa)
p

pa = get_nonsym_SNP_pa(IDs,tbl,548)
p = SNP_boxplot(dist_mat,pa)
p

pa = get_nonsym_SNP_pa_variant(IDs,tbl,'S84F')
p = SNP_boxplot(dist_mat,pa)
p

pa = get_nonsym_SNP_pa_variant(IDs,tbl,'S84F')
p = SNP_histogram(dist_mat,pa)
p

pa = get_nonsym_SNP_pa_variant(IDs,tbl,'V820A')
p = SNP_histogram(dist_mat,pa)
p
p = SNP_boxplot(dist_mat,pa)
p

p = SNP_index(dist_mat,pa)
p

test = read.table("ATCC12228_cds_SNP_sub_tables/NZ_CP022247.1_cds_WP_088922647.1_144__130372..130551.txt",sep = "\t",header=TRUE,row.names = 1,check.names=F)

sub_dist = dist(t(test))

pca_sub = prcomp(test)
plot_df = data.frame(pca_sub$x)
color_variable_factor = factor(gene_presence_df[grep('NZ_CP022247.1_cds_1804__complement_1870445',rownames(gene_presence_df)),])
color_vector = setup_colors(levels(color_variable_factor),colors)

labels = c(paste0("PC1 [",sprintf("%.1f",explvar(pca)[1]),"%]"),paste0("PC2 [",sprintf("%.1f",explvar(pca)[2]),"%]"))

p = ggplot(as.data.frame(plot_df),aes(x=PC1,y=PC2,color = color_variable_factor)) + labs(color = "gene presence") + geom_point(size=1.5, alpha=1)+ stat_ellipse(level=0.75) + scale_colour_manual(values = color_vector) + xlab(labels[1]) + ylab(labels[2]) + theme_bw()
p
