#Prerequisite= Install bedtools on your system.

#Loading needed libraries
require(data.table)
require(ggplot2)
require(dplyr)
#Loading needed functions
eval(parse("functions.R"))




###############################
#############Finding TTTT motifs
###############################
##Getting TTTT motif in both strands using seqkit
#Allow your mac to open seqkit in "System Preferences" > "Security & Privacy"
system("mkdir output") #making output folder 
if(TRUE){system("cat input/chr19.fasta | ./seqkit locate -p TTTT -m 0 -i --bed  > output/TTTT.bed")} #Put False if you have problem with seqkit
T4=fread("output/TTTT.bed",sep="\t")
T4$V2=as.integer(T4$V2+1);T4$V3=as.integer(T4$V3) #Adding 1 to V2, coz if my TTTT was from 34 to 37, then seqkit would give me from 33 to 37.



##############Choosing only those T4 that are nearby our genes(gene_info) 
######################coz there are simply too many T4s
gene_info=fread("input/genes_info.bed",sep="\t",header = T)   #Your desired set of genes
gene_info=gene_info %>% filter(V1=="chr19" & V2>40000000)
T4n=removeintersecting(T4,gene_info %>% transmute(V1,V2,V3,V4=0,V5=0,Strand),opposite = T,extra_flags = "-s")



##############Finding 3 prime signal from bedfile
if(FALSE){system("bedtools genomecov -dz -3 -ibam NascentRNASeq.bam   > 3prime_signal.bed")}
#Here we are skipping this process, since we already have a trimmed output file only for end regions of
#chromosome 19 i.e., "input/chr19_3prime_signal.bed"




############### Only taking 3prime signal nearby our chosen T4
#################coz there are simply too many 3prime signal almost billion
b=fread("input/chr19_3prime_signal.bed")    #3prime signal from the bedtools genomecov
colnames(b)=c("V1","V2","score")
b=b %>% transmute(V1,V2,V3=V2,score)
bn=removeintersecting(b,T4n %>% mutate(V2=as.integer(V2-50),V3=as.integer(V3+50)),opposite = T)



###############################
####### NOW SCORING................
###############################
a=T4n
for (vic in c(3,20,30)) #finding observed and lambda signals
{print(vic)
  a=a %>% mutate(V2=as.integer(V2-vic),V3=as.integer(V3+vic))
  dell=whatintersect(a,bn)
  dell2=dell %>% group_by_at(c(1:ncol(a))) %>% summarize(!!paste0("L",vic):=sum(score_))
  a=left_join(a,dell2)
  a[is.na(a)] = 0}
a=a %>% transmute(V1,V2=as.integer(V2+3+20+30),V3=as.integer(V3-3-20-30),V4,V5,V6,
                  observed=L3,lambda1=10*(L20-L3)/40,lambda2=10*(L30-L3)/100) #finding expected signal in 10nt as per the signal density in the neighbouting region
a=a %>% mutate(lambda=pmax(lambda1,lambda2))
a=a %>% mutate(pval=ppois(observed+10,lambda+10,lower.tail=F),adj.pval=p.adjust(pval,method="BH"))
fwrite(a,"output/T4scores.bed",col.names = T,row.names = F,sep="\t",append=F)





###############################
####### Now assigning score to our set of genes
###############################
fwrite(gene_info,"a.bed",col.names = F,row.names = F,quote = F,sep="\t",append=F)
fwrite(a %>% mutate(extra_column_to_avoid_bed12_format=0),"b.bed",col.names = F,row.names = F,quote = F,sep="\t",append=F)
system("bedtools intersect -s -loj -a a.bed -b b.bed > c.bed")
c=fread("c.bed",sep="\t")
colnames(c)=c(colnames(gene_info),paste0(colnames(a),"_t"),"delete")
u2= c %>% select(-c(delete))
u2$pval_t[u2$pval_t=="."]=1;u2$adj.pval_t[u2$adj.pval_t=="."]=1 # assigning pval of 1 to those genes which. dont intersect any T4
u2=u2 %>% mutate(T4score=-log10(as.numeric(adj.pval_t)))
u3=u2 %>% group_by(across(V1:TOTAL_INFO)) %>% slice_max(order_by = T4score, n = 1) %>% ungroup()
u3=u3 %>% group_by(across(V1:TOTAL_INFO)) %>% slice_max(order_by = as.numeric(observed_t)/as.numeric(lambda_t), n = 1) %>% 
  slice_sample(n = 1) %>% ungroup()#coz there were multiple entries
u3=u3 %>% arrange(desc(T4score))
sort(unique(u3$T4score),decreasing = T)[1:10]
u3$T4score[is.infinite(u3$T4score)]=310
sort(unique(u3$T4score),decreasing = T)[1:10]
system("rm a.bed b.bed c.bed")
fwrite(u3,"output/gene_info_with_T4score.bed",col.names = T,row.names = F,sep="\t",append=F)


#################Analyzing########
###############plotting RNA polymerase III occupancy score against T4score of snAR genes
ggplot(data=u3 %>% filter(UU=="snAR"),aes(x=pol3m,y=T4score))+
  geom_point()+theme_classic()+xlab("RNA Polymerase III occupancy score")

ggplot(data=u3 %>% filter(pol3m>2),aes(x=UU,y=T4score))+geom_boxplot(fill="#EEEEEEEE")+geom_jitter()+theme_classic()



#############################################THE END





######only extracting chr19 from hg38.fa
#awk '/^>chr19$/{p=1; print; next} p && /^>/{p=0} p' /Users/rkc/igv/genomes/seq/hg38.fa > chr19.fasta
