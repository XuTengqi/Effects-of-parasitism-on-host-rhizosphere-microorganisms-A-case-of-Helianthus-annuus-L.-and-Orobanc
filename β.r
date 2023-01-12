#########首先计算 Beta 多样性距离测度

#读取 OTU 丰度表
otu <- read.delim('otu_table.txt', row.names = 1)
otu <- t(otu)

#计算与 Beta 多样性有关的群落相异指数，例如使用 vegan 包计算 Bray-curtis 距离，详情加载 vegan 包后 ?vegdist
dis <- vegan::vegdist(otu, method = 'bray')

#以矩阵形式输出
dis <- as.matrix(dis)
write.table(dis, 'Bray-curtis.txt', sep = '\t', col.names = NA, quote = FALSE)


#########使用距离矩阵进行计算和作图

#读取 Bray-curtis 距离矩阵
dis <- read.delim('Bray-curtis.txt', row.names = 1)

#读取样本分组信息
group <- read.delim('group.txt', stringsAsFactors = FALSE)

##例如，比较 Env1、Env2、Env3 三组之间，群落的 Beta 多样性差异
#根据分组获得组内距离矩阵
env1 <- subset(group, group1 == 'Env1')$samples
dis_env1 <- dis[env1,env1]

env2 <- subset(group, group1 == 'Env2')$samples
dis_env2 <- dis[env2,env2]

env3 <- subset(group, group1 == 'Env3')$samples
dis_env3 <- dis[env3,env3]

#将矩阵转化为向量，以便用于作图和统计
dis_env1 <- as.vector(as.dist(dis_env1))
dis_env2 <- as.vector(as.dist(dis_env2))
dis_env3 <- as.vector(as.dist(dis_env3))

#构建作图数据集
dat <- data.frame(
    dis = c(dis_env1, dis_env2, dis_env3),
    group = factor(c(
        rep('Env1', length(dis_env1)), 
        rep('Env2', length(dis_env2)), 
        rep('Env3', length(dis_env3))
    ), levels = c('Env1', 'Env2', 'Env3'))
)

#使用 ggplot2 绘制各组内 Bray-curtis 距离指数分布的箱线图
library(ggplot2)

p <- ggplot(dat, aes(group, dis)) +
geom_boxplot(aes(fill = group), width = 0.6) +
scale_fill_manual(values = c('#CD5B45', '#228B22', '#00688B')) +
theme(panel.grid = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = 'black'), legend.position = 'none') +
labs(x = NULL, y = 'Bray-Curtis dissimilarity\n')

p

#三组的整体差异分析，使用 Kruskal-Wallis Test 执行，详情 ?kruskal.test
kruskal.test(dis~group, data = dat)

#如果整体显著再进行两两分组的比较，使用 Wilcoxon 秩和检验执行双侧检验，详情 ?wilcox.test
wilcox.test(dis_env1, dis_env2, alternative = 'two.sided')
wilcox.test(dis_env2, dis_env3, alternative = 'two.sided')
wilcox.test(dis_env1, dis_env3, alternative = 'two.sided')

#考虑到 Wilcoxon 秩和检验体现了中位数的差异，因此计算三组数据的中位数以评估 Beta 多样性的高低水平
median(dis_env1)
median(dis_env2)
median(dis_env3)

#基于上述统计结果，判断好组间差异后，将差异分析结果添加到箱线图中
p +
annotate('text', label = 'Kruskal-Wallis Test', x = 1, y = 0.56, size = 3) +
annotate('text', label = sprintf('italic(P) < %.3f', 0.001), x = 1, y = 0.53, size = 3, parse = TRUE) +
annotate('text', label = 'c', x = 1, y = max(dis_env1)+0.05, size = 3) +
annotate('text', label = 'a', x = 2, y = max(dis_env2)+0.05, size = 3) +
annotate('text', label = 'b', x = 3, y = max(dis_env3)+0.05, size = 3)

##例如，比较 Env1、Env2、Env3 三组内，不同季节之间群落的 Beta 多样性差异
#根据分组获得组内距离矩阵
env1_winter <- subset(group, group2 == 'Env1 winter')$samples
dis_env1_winter <- dis[env1_winter,env1_winter]
env1_summer <- subset(group, group2 == 'Env1 summer')$samples
dis_env1_summer <- dis[env1_summer,env1_summer]

env2_winter <- subset(group, group2 == 'Env2 winter')$samples
dis_env2_winter <- dis[env2_winter,env2_winter]
env2_summer <- subset(group, group2 == 'Env2 summer')$samples
dis_env2_summer <- dis[env2_summer,env2_summer]

env3_winter <- subset(group, group2 == 'Env3 winter')$samples
dis_env3_winter <- dis[env3_winter,env3_winter]
env3_summer <- subset(group, group2 == 'Env3 summer')$samples
dis_env3_summer <- dis[env3_summer,env3_summer]

#将矩阵转化为向量，以便用于作图和统计
dis_env1_winter <- as.vector(as.dist(dis_env1_winter))
dis_env1_summer <- as.vector(as.dist(dis_env1_summer))
dis_env2_winter <- as.vector(as.dist(dis_env2_winter))
dis_env2_summer <- as.vector(as.dist(dis_env2_summer))
dis_env3_winter <- as.vector(as.dist(dis_env3_winter))
dis_env3_summer <- as.vector(as.dist(dis_env3_summer))

#构建作图数据集
dat <- data.frame(
    dis = c(dis_env1_winter, dis_env1_summer, dis_env2_winter, dis_env2_summer, dis_env3_winter, dis_env3_summer), 
    group1 = factor(c(
        rep('Env1', length(c(dis_env1_winter, dis_env1_summer))), 
        rep('Env2', length(c(dis_env2_winter, dis_env2_summer))), 
        rep('Env3', length(c(dis_env3_winter, dis_env3_summer)))
    ), levels = c('Env1', 'Env2', 'Env3')), 
    group2 = factor(c(
        rep('Winter', length(dis_env1_winter)), 
        rep('Summer', length(dis_env1_summer)), 
        rep('Winter', length(dis_env2_winter)), 
        rep('Summer', length(dis_env2_summer)), 
        rep('Winter', length(dis_env3_winter)), 
        rep('Summer', length(dis_env3_summer))
        ), levels = c('Winter', 'Summer'))
)

#ggplot2 箱线图
library(ggplot2)

p <- ggplot(dat, aes(group1, dis)) +
geom_boxplot(aes(fill = group2), width = 0.5, position = position_dodge(width = 0.6)) +
scale_fill_manual(values = c('#00688B', '#CD5B45')) +
theme(panel.grid = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = 'black'), legend.key = element_blank()) +
labs(x = NULL, y = 'Bray-Curtis dissimilarity\n', fill = NULL)

p

#Env1、Env2、Env3 三组内，冬季和夏季两分组比较，使用 Wilcoxon 秩和检验执行双侧检验，详情 ?wilcox.test
wilcox.test(dis_env1_winter, dis_env1_summer, alternative = 'two.sided')
wilcox.test(dis_env2_winter, dis_env2_summer, alternative = 'two.sided')
wilcox.test(dis_env3_winter, dis_env3_summer, alternative = 'two.sided')

#计算各组数据的中位数
#由于 Env2 和 Env3 组内冬季和夏季群落的 Beta 多样性之间无差异，所以这里只计算有差异的 Env1 组
median(dis_env1_winter)
median(dis_env1_summer)

#基于上述统计结果，判断好组间差异后，将差异分析结果添加到箱线图中
p +
annotate('text', label = '———', x = 1, y = max(c(dis_env1_winter, dis_env1_summer))+0.03, size = 3) +
annotate('text', label = '———', x = 2, y = max(c(dis_env2_winter, dis_env2_summer))+0.03, size = 3) +
annotate('text', label = '———', x = 3, y = max(c(dis_env3_winter, dis_env3_summer))+0.03, size = 3) +
annotate('text', label = '***', x = 1, y = max(c(dis_env1_winter, dis_env1_summer))+0.05, size = 3) +
annotate('text', label = 'none', x = 2, y = max(c(dis_env2_winter, dis_env2_summer))+0.05, size = 3) +
annotate('text', label = 'none', x = 3, y = max(c(dis_env3_winter, dis_env3_summer))+0.05, size = 3)
