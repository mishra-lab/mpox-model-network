# some plotting tools

suppressPackageStartupMessages({
  library('ggplot2')
  library('viridis')
})

.plot.ext = 'pdf'

fig.save = function(fname,quoted=NULL,g=last_plot(),w=7,h=7,margin=c(0,0,0,0)){
  # usually ggsave, but alternatively save base figure like fig.save('a',quote(plot(1:3)))
  fname = paste0(fname,'.',.plot.ext)
  print(paste('saving:',fname))
  if (is.null(quoted)){ # usually ggplot with last_plot()
    ggsave(fname,plot=g,w=w,h=h)
  } else { # using grDevices wrapping of quoted
    if (.plot.ext != 'pdf'){ w=w*72; h=h*72 }
    get(.plot.ext)(fname,width=w,height=h)
    par(mai=c(h,w,h,w)*margin) # adjust margins
    eval(quoted) # build the plot
    dev.off()
  }
}
