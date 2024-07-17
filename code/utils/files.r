# some tools for file stuff

root.path = function(...,create=FALSE){
  # make a file path starting one level above /code/
  # e.g. root.path('foo','bar') -> '.../{projectroot}/foo/bar'
  root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create & !dir.exists(dirname(path))){ dir.create(dirname(path),recursive=TRUE) }
  return(path)
}

tmp.out = function(out,slug,path='.tmp',ext='.out'){
  # print some text (out) and also save it to a temporary file; also, path can be a vector
  # e.g. path=c('a','b') -> 'a/b/{slug}-1970-01-01.out'
  file.out = kw.call(file.path,c(as.list(path),paste0(slug,'-',Sys.Date(),ext)))
  print(out); sink(file.out); print(out); sink()
}
