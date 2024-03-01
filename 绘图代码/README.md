# 绘图代码
## 1. 添加抖动的柱状图（蜂群柱状图）
需要加载的包：
````R
library(tidyverse)
library(ggbeeswarm)
library(ggsignif)
library(plyr)
````
需要准备的数据：数据框
代码中为了示例已经生成一个300行2列的长数据dt_long和3行、3列的data.fram:dt_sum:

示例图片：  