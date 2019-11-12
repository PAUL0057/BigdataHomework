library(ggplot2)
data = read.csv('d:\\大数据相关\\抽样\\LoanStats3c\\LoanStats3c.csv')

income = na.omit(data[,'annual_inc'])

##数据可视化观察数据情况
plot(sort(income))
summary(income)
#数据严重右偏
my_break = c(0,10000*(1:10),1500000,100000*(2:5),max(income))
income_cut = cut(income,breaks=my_break)
table(income_cut)
freq = table(income_cut)/length(income)
mydata = data.frame(table(income_cut))
ggplot(data = mydata , aes(x = income_cut,y = Freq,fill = income_cut))+
  geom_bar(stat = 'identity')+
  labs(title="income条形图",x='income',y='count')+
  theme(axis.text.x = element_text(angle = 295, hjust = 0.5, vjust = 0.5)) 

#通过指数分布选取68个样本容量
x = seq(10,20,0.1)
y = 2^(x)
sampl = round(y)[1:68]
sample1 = function(x){
  y = sample(income,x)
  y = c(y,rep(NA,sampl[length(sampl)]-length(y)))
  return(y)
}
sample_num = as.matrix(sampl)
sample_result = apply(sample_num,1,sample1)


#计算样本质量过程
get_score = function(sample_result){
  income_sample = cut(na.omit(sample_result),breaks = my_break)
  freq1 = table(income_sample)/length(na.omit(sample_result))+1e-15#防止生成Inf
  score = exp(-sum((freq1-freq)*(log(freq1/freq))))
  return(score)
}

#随机抽样
scores = NULL
for ( i in 1:50) {
  sample_result = apply(sample_num,1,sample1)
  score = apply(sample_result,2,get_score)
  scores = rbind(scores,score)
}
scores_final = apply(scores,2,mean)



#分层抽样
sample_num_table = sampl %*% t(as.vector(freq))
sample_data = data.frame(income,income_cut)
box = levels(income_cut)
sample2 = function(i,j){# i 是68选1 j是16个分层

  x = sample_num_table[i,j]
  box1 = box[j]
  y = sample(sample_data[sample_data[,2]==box1,]$income,x)

  return(y)
}

sample_result_func = function(){
z = NULL
for (i in 1:68) {
  y = NULL
  for (j in 1:16) {
    y1 = sample2(i,j)
    y = c(y,y1)
  }
  y = c(y,rep(NA,(sampl[68]-length(y))))
  z = cbind(z,y)
}
return (z)
}

scores = NULL
for ( i in 1:50) {
  sample_result = sample_result_func()
  score = apply(sample_result,2,get_score)
  scores = rbind(scores,score)
}
scores_final = apply(scores,2,mean)
