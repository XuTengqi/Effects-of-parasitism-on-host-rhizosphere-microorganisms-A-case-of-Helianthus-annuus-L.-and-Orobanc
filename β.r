#########���ȼ��� Beta �����Ծ�����

#��ȡ OTU ��ȱ�
otu <- read.delim('otu_table.txt', row.names = 1)
otu <- t(otu)

#������ Beta �������йص�Ⱥ������ָ��������ʹ�� vegan ������ Bray-curtis ���룬������� vegan ���� ?vegdist
dis <- vegan::vegdist(otu, method = 'bray')

#�Ծ�����ʽ���
dis <- as.matrix(dis)
write.table(dis, 'Bray-curtis.txt', sep = '\t', col.names = NA, quote = FALSE)


#########ʹ�þ��������м������ͼ

#��ȡ Bray-curtis �������
dis <- read.delim('Bray-curtis.txt', row.names = 1)

#��ȡ����������Ϣ
group <- read.delim('group.txt', stringsAsFactors = FALSE)

##���磬�Ƚ� Env1��Env2��Env3 ����֮�䣬Ⱥ��� Beta �����Բ���
#���ݷ��������ھ������
env1 <- subset(group, group1 == 'Env1')$samples
dis_env1 <- dis[env1,env1]

env2 <- subset(group, group1 == 'Env2')$samples
dis_env2 <- dis[env2,env2]

env3 <- subset(group, group1 == 'Env3')$samples
dis_env3 <- dis[env3,env3]

#������ת��Ϊ�������Ա�������ͼ��ͳ��
dis_env1 <- as.vector(as.dist(dis_env1))
dis_env2 <- as.vector(as.dist(dis_env2))
dis_env3 <- as.vector(as.dist(dis_env3))

#������ͼ���ݼ�
dat <- data.frame(
    dis = c(dis_env1, dis_env2, dis_env3),
    group = factor(c(
        rep('Env1', length(dis_env1)), 
        rep('Env2', length(dis_env2)), 
        rep('Env3', length(dis_env3))
    ), levels = c('Env1', 'Env2', 'Env3'))
)

#ʹ�� ggplot2 ���Ƹ����� Bray-curtis ����ָ���ֲ�������ͼ
library(ggplot2)

p <- ggplot(dat, aes(group, dis)) +
geom_boxplot(aes(fill = group), width = 0.6) +
scale_fill_manual(values = c('#CD5B45', '#228B22', '#00688B')) +
theme(panel.grid = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = 'black'), legend.position = 'none') +
labs(x = NULL, y = 'Bray-Curtis dissimilarity\n')

p

#�����������������ʹ�� Kruskal-Wallis Test ִ�У����� ?kruskal.test
kruskal.test(dis~group, data = dat)

#������������ٽ�����������ıȽϣ�ʹ�� Wilcoxon �Ⱥͼ���ִ��˫����飬���� ?wilcox.test
wilcox.test(dis_env1, dis_env2, alternative = 'two.sided')
wilcox.test(dis_env2, dis_env3, alternative = 'two.sided')
wilcox.test(dis_env1, dis_env3, alternative = 'two.sided')

#���ǵ� Wilcoxon �Ⱥͼ�����������λ���Ĳ��죬��˼����������ݵ���λ�������� Beta �����Եĸߵ�ˮƽ
median(dis_env1)
median(dis_env2)
median(dis_env3)

#��������ͳ�ƽ�����жϺ�������󣬽�������������ӵ�����ͼ��
p +
annotate('text', label = 'Kruskal-Wallis Test', x = 1, y = 0.56, size = 3) +
annotate('text', label = sprintf('italic(P) < %.3f', 0.001), x = 1, y = 0.53, size = 3, parse = TRUE) +
annotate('text', label = 'c', x = 1, y = max(dis_env1)+0.05, size = 3) +
annotate('text', label = 'a', x = 2, y = max(dis_env2)+0.05, size = 3) +
annotate('text', label = 'b', x = 3, y = max(dis_env3)+0.05, size = 3)

##���磬�Ƚ� Env1��Env2��Env3 �����ڣ���ͬ����֮��Ⱥ��� Beta �����Բ���
#���ݷ��������ھ������
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

#������ת��Ϊ�������Ա�������ͼ��ͳ��
dis_env1_winter <- as.vector(as.dist(dis_env1_winter))
dis_env1_summer <- as.vector(as.dist(dis_env1_summer))
dis_env2_winter <- as.vector(as.dist(dis_env2_winter))
dis_env2_summer <- as.vector(as.dist(dis_env2_summer))
dis_env3_winter <- as.vector(as.dist(dis_env3_winter))
dis_env3_summer <- as.vector(as.dist(dis_env3_summer))

#������ͼ���ݼ�
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

#ggplot2 ����ͼ
library(ggplot2)

p <- ggplot(dat, aes(group1, dis)) +
geom_boxplot(aes(fill = group2), width = 0.5, position = position_dodge(width = 0.6)) +
scale_fill_manual(values = c('#00688B', '#CD5B45')) +
theme(panel.grid = element_blank(), panel.background = element_blank(), 
    axis.line = element_line(colour = 'black'), legend.key = element_blank()) +
labs(x = NULL, y = 'Bray-Curtis dissimilarity\n', fill = NULL)

p

#Env1��Env2��Env3 �����ڣ��������ļ�������Ƚϣ�ʹ�� Wilcoxon �Ⱥͼ���ִ��˫����飬���� ?wilcox.test
wilcox.test(dis_env1_winter, dis_env1_summer, alternative = 'two.sided')
wilcox.test(dis_env2_winter, dis_env2_summer, alternative = 'two.sided')
wilcox.test(dis_env3_winter, dis_env3_summer, alternative = 'two.sided')

#����������ݵ���λ��
#���� Env2 �� Env3 ���ڶ������ļ�Ⱥ��� Beta ������֮���޲��죬��������ֻ�����в���� Env1 ��
median(dis_env1_winter)
median(dis_env1_summer)

#��������ͳ�ƽ�����жϺ�������󣬽�������������ӵ�����ͼ��
p +
annotate('text', label = '������', x = 1, y = max(c(dis_env1_winter, dis_env1_summer))+0.03, size = 3) +
annotate('text', label = '������', x = 2, y = max(c(dis_env2_winter, dis_env2_summer))+0.03, size = 3) +
annotate('text', label = '������', x = 3, y = max(c(dis_env3_winter, dis_env3_summer))+0.03, size = 3) +
annotate('text', label = '***', x = 1, y = max(c(dis_env1_winter, dis_env1_summer))+0.05, size = 3) +
annotate('text', label = 'none', x = 2, y = max(c(dis_env2_winter, dis_env2_summer))+0.05, size = 3) +
annotate('text', label = 'none', x = 3, y = max(c(dis_env3_winter, dis_env3_summer))+0.05, size = 3)
