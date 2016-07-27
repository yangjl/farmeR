if(!require(QuickShinyApp)) devtools::install_github('deruncie/QuickShiny')


slider_params = list(
  d = c(0,-1.5,1.5),
  q = c(0.5,0.01,0.99)
)

plot_fun = function(d,q) {
  area = 0.25^2*pi
  p = 1-q
  a = 1
  alpha = a + d*(q-p)
  mean = a*(p-q) + 2*d*p*q
  plot(NA,NA,xlim = c(-0.5,2.5),ylim = c(-a,a)*1.5,xlab = '# A1',ylab = 'genetic value')
  symbols(c(0,1,2),y=c(-1,d,1),circles = sqrt(area*c(q^2,2*p*q,p^2)/pi),add=T,inches=F,bg = 'grey70')
  abline(mean-2*p*alpha,alpha)
  text(1,d - 0.5*c(1,-1)[(d<0)+1],sprintf('alpha = %0.2f',alpha))
}

call = 'plot_fun(d,q)'
library(QuickShinyApp)
run_shiny(call,slider_params)
