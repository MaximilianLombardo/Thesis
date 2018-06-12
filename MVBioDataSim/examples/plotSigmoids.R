library(ggplot2)
library(grid)
library(gridExtra)

#fix theta, variate mu
theta=0.3
mu=c(2, 6, 10)

x=seq(from = 0, to = 1.2, length.out = 100)

mmu=mu[2]
mumin=mu[1]
mumax=mu[3]
ymu=x^mmu/(x^mmu+theta^mmu)
ymmin=x^mumin/(x^mumin+theta^mumin)
ymmax=x^mumax/(x^mumax+theta^mumax)

#fix h, variate theta
theta=c(0.3, 0.5, 0.8)
mu=5

tt=theta[2]
tmin=theta[1]
tmax=theta[3]
yt=x^mu/(x^mu+tt^mu)
ytmin=x^mu/(x^mu+tmin^mu)
ytmax=x^mu/(x^mu+tmax^mu)
df=data.frame(x=x, ymu=ymu, ymmax=ymmax, ymmin=ymmin, yt=yt, ytmin=ytmin, ytmax=ytmax)

alpha=0.2

mfixed=ggplot(df, aes(x=x))+geom_ribbon(aes(y=ymu, ymin=ymmin, ymax=ymmax),alpha=alpha, fill='blue') + 
  geom_line(aes(y=ymu), colour='blue') + ylab(expression(h(x ~ ";" ~ list(theta, h))))

tfixed=ggplot(df, aes(x=x))+geom_ribbon(aes(y=yt, ymin=ytmin, ymax=ytmax), alpha=alpha, fill='red')+
  geom_line(aes(y=yt),colour='red') + ylab(expression(h(x ~ ";" ~ list(theta, h))))

grob=arrangeGrob(tfixed, mfixed, ncol=2)
print(grob)
#ggsave(filename = 'hill.pdf', plot = grob)
