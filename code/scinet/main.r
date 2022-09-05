# load
source('model/do/vt0/main.r')
# config
.n.cores = 80
vt0$obj1$N.s  = 100
vt0$obj2a$N.s = 100
vt0$obj2b$N.s = 100
# run
vt0.obj('obj1', .run=TRUE)
vt0.obj('obj2a',.run=TRUE)
vt0.obj('obj2b',.run=TRUE)
