


print("The functions are :")



#######################################
#######################################
#################### REMOVE INTERSECTING GENES
#######################################
print("removeintersecting")
removeintersecting=function(a,b,opposite=F,extra_flags="")
{
  flags=ifelse(opposite==F,"-v","-u")
  print("Writing files to perform bedtools intersect.......")
  write.table(a ,"a.bed",sep="\t",row.names = F,col.names = F,quote=F,append=F)
  write.table(b ,"b.bed",sep="\t",row.names = F,col.names = F,quote=F,append=F)
  flags=paste(flags,extra_flags)
  code=paste("bedtools intersect",flags,"-a a.bed -b b.bed > c.bed")
  print(paste("Running -->",code))
  system(code)
  print("Now reading the intersected file and shortly returning it.....")
  c=fread("c.bed",sep="\t")
  colnames(c)=colnames(a)
  system("rm a.bed b.bed c.bed")
  return(c)
}


#######################################
#######################################
#################### What genes intersect
#######################################
print("whatintersect")
whatintersect=function(a,b,cola=NULL,colb=NULL,f=NULL)
{
  fwrite(a,"a.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)
  fwrite(b,"b.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)
  if(is.null(f)) {system("bedtools intersect -wa -wb -a a.bed -b b.bed > c.bed ")} else
  {system(paste0("bedtools intersect -wa -wb -f ",f," -a a.bed -b b.bed > c.bed "))}
  c=fread("c.bed",sep="\t")
  if(is.null(cola)){cola=colnames(a)}
  if(is.null(colb)){colb=paste0(colnames(b),"_")}
  colnames(c)=c(cola,colb)
  system("rm a.bed b.bed c.bed")
  return(c)
}





#######################################
#######################################
#################### bedtools
#######################################
print('bedtools(a,b,commands="merge -flag1 -flag2 ..... -a a.bed -b bed ")')
bedtools=function(a=NULL,b=NULL,commands=NULL){  if(is.null(commands)) {print("No commands and flags were given"); break}
  print("Writing files in progress ....")  if(!is.null(a)){fwrite(a,"a.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)}  if(!is.null(b)){fwrite(b,"b.bed",row.names = F,col.names = F,quote=F,sep="\t",append=F)}  code=paste0("bedtools ",commands," > c.bed ")
  print("Performing")
  print(code)  system(code)
  print("Now reading and returning the output file ...")  c=fread("c.bed",sep="\t")  system("rm a.bed b.bed c.bed")  return(c)}

