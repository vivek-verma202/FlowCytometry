logicleTrans=function(Frame,Markers){
x=exprs(Frame)
x2=x[,Markers]
lgcl=logicleTransform(w = 0.5,t=262143,m=4.5)
trans=transformList(colnames(exprs(Frame)),lgcl)
after=transform(Frame,trans)
x[,Markers]=matrix(exprs(after)[,Markers],dim(x2))
exprs(Frame)=x
return(Frame)
}
