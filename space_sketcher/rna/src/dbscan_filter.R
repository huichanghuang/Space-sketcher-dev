##This script is used to perform dbscan clustering for each cell barcode.
##And then calculate the ratio of single-cluster, multi-cluster and no-cluster.
##In single-cluster, calculate the weighted centroid of the cluster.
##In multi-cluster, choose the max cluster as the cluster, and calculate the weighted centroid of the cluster.
##Filter out no-cluster cell barcodes.
library(ggplot2)
library(dplyr)
library(readr)
library(dbscan)
library(gridExtra)
library(argparse)


# Define the command line arguments
parser <- ArgumentParser()
parser$add_argument("--infile", type = "character", help = "The directory of the matrix")
parser$add_argument("--samplename", type = "character", help = "The name of the sample")
parser$add_argument("--outdir", type = "character", help = "The output directory")
parser$add_argument("--maxumi", type = "numeric", help = "The number of max umi count of spatial barcode to filter out, default: 5000", default = 5000)
parser$add_argument("--minumi", type = "numeric", help = "The number of min umi count of spatial barcode to filter out, default: 2", default = 2)
args <- parser$parse_args()

infile <- args$infile
samplename <- args$samplename
outdir <- args$outdir
maxumi <- args$maxumi
minumi <- args$minumi

# rm(list=ls())
# infile <- "/data01/huanghuichang/spatialRNA_analysis/20250120/N16-puck-9-test/N16-puck-9.CB_coord.total.txt"
# samplename <- "N16-puck-9-filterumi-lowandhigh"
# outdir <- "/data01/huanghuichang/spatialRNA_analysis/20250120/N16-puck-9-test"

#Load the data
CB_coord_total <- read_delim(infile, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
head(CB_coord_total)

###group umi count by subsb
subsb_umi_summary <- CB_coord_total %>%
  group_by(subsb) %>%
  summarise(umi_count_sum = sum(umi_count))
###sort
subsb_umi_summary <- subsb_umi_summary %>%
  arrange(desc(umi_count_sum))
###add rank
subsb_umi_summary <- subsb_umi_summary %>%
  mutate(rank = row_number())
head(subsb_umi_summary)

total_valid_sb <- dim(subsb_umi_summary)[1]
sb_umi_count_larger_5000 <- dim(subsb_umi_summary[subsb_umi_summary$umi_count_sum>5000,])[1]
sb_umi_count_smaller_2 <- dim(subsb_umi_summary[subsb_umi_summary$umi_count_sum<2,])[1]

##spatial barcode - umi 骤降图
p <- ggplot(subsb_umi_summary, aes(x = rank, y = umi_count_sum)) +
  geom_point() +  # 绘制散点图
  geom_line() +   # 绘制连接线
  geom_hline(yintercept = 5000, linetype = "dashed", color = "red") +  
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +  
  scale_x_log10() +  # 设置x轴为对数刻度
  scale_y_log10() +  # 设置y轴为对数刻度
  labs(
    title = "Log-Log Plot of UMI Count Sum",
    x = "Rank",
    y = "UMI Count Sum"
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(
    panel.background = element_rect(fill = "white"),  # 设置面板背景为白色
    plot.background = element_rect(fill = "white", color = "white")  # 设置整个图形背景为白色
  )

ggsave(paste0(outdir,"/", samplename, ".spatial_umi_knee_plot.png"), plot = p, width = 8, height = 6, dpi = 300)
###output subsb_umi_summary
subsb_umi_summary_out <- paste0(outdir,"/", samplename, ".spatial_umi_knee_plot.txt")
write.table(subsb_umi_summary, file = subsb_umi_summary_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


### 筛选出umi_count_sum大于5000或小于2的subsb
subsb_to_filter <- subsb_umi_summary %>%
  filter(umi_count_sum > 5000 | umi_count_sum < 2) %>%
  select(subsb)  # 仅保留subsb列

CB_coord_total_filtered <- CB_coord_total %>%
  filter(!subsb %in% subsb_to_filter$subsb)

####dbscan 聚类
###unique cell barcodes
unique_barcodes <- unique(CB_coord_total_filtered$cb)

no_cluster <- 0
single_cluster <- 0
multi_cluster <- 0
filtereddf <- data.frame()
cb_cluster <- data.frame()
###clustering for each cell barcode
for (barcode in unique_barcodes) {

  cb_sel <- CB_coord_total_filtered %>% filter(cb == barcode) %>% arrange(umi_count)
  ###check if there are enough data points to perform clustering
  if (nrow(cb_sel) < 6) {
    cb_sel$db_clu <- as.factor(0)
    no_cluster <- no_cluster + 1
    clustertype <- "no-cluster"
    cb_cluster_temp <- data.frame(cb = barcode, cluster = clustertype)
    cb_cluster <- rbind(cb_cluster, cb_cluster_temp)
    next
  }

  db <- dbscan(cb_sel[, c("xcoord", "ycoord")], eps = 150, minPts = 6, weights = sqrt(cb_sel$umi_count))
  cb_sel$db_clu <- as.factor(db$cluster)
  ######detect multi-clusters or single-cluster
  head(cb_sel)
  cluster_table <- table(cb_sel$db_clu)
  cluster_counts <- as.numeric(cluster_table)
  ###remove the first cluster count cause it is the noise cluster
  cluster_counts <- cluster_counts[-1]
  if (length(cluster_counts) > 0) {
    maxcluser <- max(cluster_counts)
    if (length(cluster_counts) > 1) {
      second_max <- sort(cluster_counts, decreasing = TRUE)[2]
    } else {
      second_max <- 0
    }
    if (second_max > maxcluser/2) {
      multi_cluster <- multi_cluster + 1
      clustertype <- "multi-cluster"
    } else {
      single_cluster <- single_cluster + 1
      clustertype <- "single-cluster"
    }
    ###Weighted centroid calculation
    cluster1_cell_cb <- cb_sel %>% filter(db_clu == which.max(cluster_counts))
    xcor <- cluster1_cell_cb$xcoord
    ycor <- cluster1_cell_cb$ycoord
    weights <- cluster1_cell_cb$umi_count
    weighted_centroid_x <- round(sum(xcor * weights) / sum(weights),0)
    weighted_centroid_y <- round(sum(ycor * weights) / sum(weights),0)
    weighted_centroid_cell <- data.frame(cb = barcode, xcoord = weighted_centroid_x, ycoord = weighted_centroid_y)
    filtereddf <- rbind(filtereddf, weighted_centroid_cell)
  } else {
    no_cluster <- no_cluster + 1
    clustertype <- "no-cluster"
  }
  cb_cluster_temp <- data.frame(cb = barcode, cluster = clustertype)
  cb_cluster <- rbind(cb_cluster, cb_cluster_temp)
}

###calculating the ratio of single-cluster, multi-cluster and no-cluster
cluster_counts <- c(single_cluster, multi_cluster, no_cluster)
cluster_ratio <- round(cluster_counts*100/sum(cluster_counts), 2)
###output statistics
# output_file <- file("barcode_dbscan_cluster_distribution.txt", "w")
output_file <- file(paste(outdir, "/", samplename, ".barcode_dbscan_cluster_distribution.txt", sep = ""), "w")
writeLines("clusternum\tcount\tratio", output_file)
writeLines(paste("single-cluster\t", cluster_counts[1], "\t", cluster_ratio[1]), output_file)
writeLines(paste("multi-cluster\t", cluster_counts[2], "\t", cluster_ratio[2]), output_file)
writeLines(paste("no-cluster\t", cluster_counts[3], "\t", cluster_ratio[3]), output_file)
writeLines(paste("total\t", sum(cluster_counts), "\t", 1), output_file)

head(filtereddf)
outfile <- paste(outdir, "/", samplename, ".CB_coord.dbscan_filtered.txt", sep = "")
write.table(filtereddf, file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###unique cell barcodes
cb_cluster <- unique(cb_cluster)
outcbcluster <- paste(outdir, "/", samplename, ".CB_cluster.txt", sep = "")
write.table(cb_cluster, file = outcbcluster, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# 将clustertype中的cluster信息合并到CB_coord_total_filtered中
CB_coord_total_filtered <- CB_coord_total_filtered %>%
  left_join(cb_cluster, by = "cb")
head(CB_coord_total_filtered)

# 按cb分组计算每个细胞中top100 sb的umi count 的均值
cb_umi_mean_top100 <- CB_coord_total_filtered %>%
  group_by(cb) %>%
  arrange(desc(umi_count)) %>%  # 按umi_count降序排序
  slice_head(n = 100) %>%  # 选择每个cb组的前100个sb
  group_by(cb, cluster) %>%  # 重新按cb分组
  summarise(umi_count_mean = mean(umi_count)) %>%  # 计算前100个sb的umi_count均值
  ungroup()  # 去掉分组

head(cb_umi_mean_top100)
# 按umi_count_mean降序排序
cb_umi_mean_top100 <- cb_umi_mean_top100 %>%
  arrange(desc(umi_count_mean))

# 添加序号列
cb_umi_mean_top100 <- cb_umi_mean_top100 %>%
  mutate(rank = row_number())
head(cb_umi_mean_top100)

median_top100_umi_cell_mean <- median(cb_umi_mean_top100$umi_count_mean)
mean_top100_umi_cell_mean <- mean(cb_umi_mean_top100$umi_count_mean)

# 绘制对数-对数曲线
p <- ggplot(cb_umi_mean_top100, aes(x = rank, y = umi_count_mean, color = cluster)) +
  geom_point() +  # 绘制散点图
  geom_line() +   # 绘制连接线
  scale_x_log10() +  # 设置x轴为对数刻度
  scale_y_log10() +  # 设置y轴为对数刻度
  labs(
    title = "Log-Log Plot of UMI Count Mean by CB",
    x = "Rank",
    y = "Mean of UMI Count of Each CB",
    color = "Cluster"
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(
    panel.background = element_rect(fill = "white"),  # 设置面板背景为白色
    plot.background = element_rect(fill = "white", color = "white")  # 设置整个图形背景为白色
  )

ggsave(paste0(outdir, "/", samplename, ".cb_umi_mean_knee_plot.png"), plot = p, width = 8, height = 6, dpi = 300)
###output cb_umi_mean_top100
cb_umi_mean_top100_out <- paste0(outdir,"/", samplename, ".cb_umi_mean_knee_plot.txt")
write.table(cb_umi_mean_top100, file = cb_umi_mean_top100_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###output new stat
stats <- c(
  samplename = samplename,
  total_valid_sb = total_valid_sb,
  sb_umi_count_larger_5000 = sb_umi_count_larger_5000,
  sb_umi_count_smaller_2 = sb_umi_count_smaller_2,
  median_top100_umi_cell_mean = median_top100_umi_cell_mean,
  mean_top100_umi_cell_mean = mean_top100_umi_cell_mean
)
stats_matrix <- matrix(stats, nrow = 1, byrow = TRUE)
colnames(stats_matrix) <- names(stats)

stats_matrix_out <- paste0(outdir,"/", samplename, ".sb_cb_umi_stats.txt")
write.table(stats_matrix, file = stats_matrix_out, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

###cell barcode vs spatial barcodes 种类数
# 将clustertype中的cluster信息合并到CB_coord_total中
CB_coord_total <- CB_coord_total %>%
  left_join(cb_cluster, by = "cb")
head(CB_coord_total)

sb_count <- CB_coord_total %>%
  group_by(cb,cluster) %>%
  summarise(sbcounts = n_distinct(sb)) %>%
  ungroup()


# 根据种类数量从大到小排序，并新增序号列
sb_count <- sb_count %>%
  arrange(desc(sbcounts)) %>%
  mutate(rank = row_number())
head(sb_count)

# 绘制对数-对数图
p <- ggplot(sb_count, aes(x = rank, y = sbcounts, color = cluster)) +
  geom_point() +  # 绘制散点图
  geom_line() +   # 绘制连接线
  scale_x_log10() +  # 设置x轴为对数刻度
  scale_y_log10() +  # 设置y轴为对数刻度
  labs(
    title = "Log-Log Plot of sb counts by CB",
    x = "Rank",
    y = "sb counts"
  ) +
  theme_minimal() +  # 使用简洁主题
  theme(
    panel.background = element_rect(fill = "white"),  # 设置面板背景为白色
    plot.background = element_rect(fill = "white", color = "white")  # 设置整个图形背景为白色
  )

ggsave(paste0(outdir, "/", samplename, ".cb_sb_counts.png"), plot = p, width = 8, height = 6, dpi = 300)
