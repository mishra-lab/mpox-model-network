# some plotting tools

suppressPackageStartupMessages({
  library('ggplot2')
  library('scales')
  library('ggExtra')
  library('ggridges')
  library('viridis')
})

.plot.ext = 'pdf'

fig.save = function(...,g=last_plot(),w=7,h=7,ext=.plot.ext){
  # wrapper for ggsave
  fname = paste0(...,'.',ext)
  print(paste('saving:',fname))
  ggsave(fname,plot=g,w=w,h=h)
}

nsqrt_trans = function(){
  trans_new('nsqrt',
    function(x){sign(x)*sqrt(abs(x))},
    function(x){x^2*sign(x)})
}
