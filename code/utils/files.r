# some tools for file stuff

file.ext = function(fname){
  # e.g. 'apple.pie.zip' -> 'zip'
  fname.split = strsplit(fname,'\\.')
  return(fname.split[[length(fname.split)]])
}

root.path = function(...,create=FALSE){
  # make a file path starting one level above /code/
  # e.g. root.path('foo','bar') -> '.../{projectroot}/foo/bar'
  root = strsplit(file.path(getwd(),''),file.path('','code',''))[[1]][1]
  path = file.path(root,...)
  if (create & !dir.exists(dirname(path))){ dir.create(dirname(path),recursive=TRUE) }
  return(path)
}