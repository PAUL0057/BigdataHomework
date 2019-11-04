library(nycflights13)
data = flights
data1 = data[data$arr_delay>=2,]
data2 = data1[data1$dest%in% c('IAH','HOU'),]
data3 = data2[data2$carrier%in%c('UA','DL'),]
data4 = data3[data3$month%in%c(7,8,9),]
data5 = data[(data$arr_delay>120)&(data$dep_delay)<=0,]
data6 = data[(data$dep_delay>=60)&(data$arr_delay<30),]
data7 = data[data$dep_time=<600 & data$dep_time>0,]

#3
delay_max = max(data$dep_delay,na.rm=T)
data_dep_delay = data[!is.na(data$dep_delay),]
data_dep_delay[data_dep_delay$dep_delay==1301,]$flight
earlest_time = min(data[data$flight==51,]$dep_time)
print(earlest_time)

#4
num_time = function(x)
  {
  if (is.na(x)==T){y = NA}
  else if (x <=60 ) {y = x
  } else  {
    x = as.character(x)
    minu = as.numeric(substr(x,nchar(x)-1,nchar(x)))
    hour = as.numeric(substr(x,1,nchar(x)-2))
    y = hour*60+minu}
  return(y)
}
data$dep_time_sec = 0
for(i in 1:length(data$dep_time))
    {data$dep_time_sec[i] = num_time(data$dep_time[i])}
data$sched_dep_time_sec = 0
for(i in 1:length(data$dep_time))
{data$sched_dep_time_sec[i] = num_time(data$sched_dep_time[i])}

library(magrittr)
dep_delay = data$dep_delay
not_cancelled <- flights %>%
  filter(!is.na(dep_delay), !is.na(arr_delay))
tapply(not_cancelled$distance,not_cancelled$dest,sum)
tapply(not_cancelled$dep_delay,not_cancelled$dest,mean)
table(not_cancelled$dest)
rm_table <- not_cancelled %>%
  filter(dest != 'HNL',origin !='HNL')


#|(6)group the data by each day, then
#|a.compute the average arriving delay time
#|b.When do the first and last flights leave each day?


tapply(not_cancelled$arr_delay,not_cancelled$day,mean)
tapply(not_cancelled$dep_time,not_cancelled$day,min)
tapply(not_cancelled$dep_time,not_cancelled$day,max)

airport =airports[,c('faa','lat','lon')]
a = left_join(not_cancelled,airport,c('origin'='faa'))
a$origin_lon = a$lon
a$origin_lat = a$lat
b = left_join(not_cancelled,airport,c('dest'='faa'))
a$dest_lon = b$lon
a$dest_lat = b$lat

plane = planes 
c = left_join(not_cancelled,planes,c('tailnum'='tailnum'))
c = na.omit(c)
cor.test(c$dep_delay,c$year.y)

head(weather)
d = left_join(not_cancelled,weather,c('year'='year','month'='month','day'='day','hour'='hour','origin'='origin'))
library(corrplot)
d = na.omit(d)
corrplot(d$dep_delay,d$temp,d$dewp)
cor(d$dep_delay,d[,c('temp','dewp','humid','wind_speed','wind_gust')])

table_tail = table(flights$tailnum)
flight_name = names(table_tail[table_tail>100])
flights[flights$tailnum %in% flight_name,]
