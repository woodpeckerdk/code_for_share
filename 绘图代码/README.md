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
<img src="https://github.com/woodpeckerdk/my_r_code/blob/main/%E7%BB%98%E5%9B%BE%E4%BB%A3%E7%A0%81/%E7%A4%BA%E4%BE%8B%E5%9B%BE%E7%89%87/bar_errorbar_point.png" width="20%">