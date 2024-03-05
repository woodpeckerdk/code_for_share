# Load package
rm(list = ls())
library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
library(plyr)
# create data
data <- data.frame(
  Ctrl = abs(rnorm(100,4.5,0.3)),
  group_1 = abs(rnorm(100,4,0.3)),
  group_2 =abs(rnorm(100,4.2,0.5))
)
data <- data[1:10,]
write.csv(data,"分组数据示例.csv")
data<- read.csv("分组数据示例.csv")
head(data)
# wide to long
dt_long <- gather(data,group,value)#变换为长数据
head(dt_long)
str(dt_long)
dt_long$group <- factor(dt_long$group,levels = c('Ctrl','group_1','group_2'))
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable to be summariezed
# groupnames : vector of column names to be used as grouping variables
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
dt_sum <- data_summary(dt_long, varname="value",groupnames="group")
class(dt_sum)
# draw
ggplot() + 
  geom_errorbar(dt_sum,mapping = aes(x=group, y=value,ymin=value-sd, ymax=value+sd), 
                width=.4,position=position_dodge(.8)) +
  geom_signif(dt_long,mapping = aes(x=group, y=value),
              comparisons = list(c("Ctrl","CD"),c("Ctrl","UC")),
              test = "t.test",
              step_increase = 0.2,
              tip_length = 0,
              textsize = 6,
              size = 1,
              map_signif_level = F)+
  geom_beeswarm(dt_long,mapping = aes(x=group, y=value, fill=group),
              shape = 21,color = 'black',size = 5,cex = 1.2,stroke = 0.6)+
  scale_fill_manual(values=c('#73bbaf','#d15e67','#6c43a6'))+
  scale_y_continuous(expand = c(0,0),limits = c(3,5.5))+
  labs(x = "", y = "Shannon index")+
  theme_classic()+
  theme(axis.line.x=element_line(color="black",size=0.8),
        axis.line.y=element_line(color="black",size=0.8),
        axis.ticks.x=element_line(color="black",size=0.8),
        axis.ticks.y=element_line(color="black",size=0.8),
        axis.text.x = element_text(color="black",size=14),
        axis.title.y=element_text(color="black",size=14))
ggsave('bar_errorbar_point.png',width = 5,height = 6)
