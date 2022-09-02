cat $1.tex > which.aux
pdflatex main ;\
pdfcrop --margins '5 5 5 5' main.pdf ../../out/fig/tikz/$1.pdf