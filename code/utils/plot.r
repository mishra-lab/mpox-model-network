# some plotting tools

suppressPackageStartupMessages({
  library('ggplot2')
  library('viridis')
})

.plot.ext = 'pdf'

fig.save = function(...,g=last_plot(),w=7,h=7){
  # wrapper for ggsave
  fname = paste0(...,'.',.plot.ext)
  print(paste('saving:',fname))
  ggsave(fname,plot=g,w=w,h=h)
}
