library(ape)
  
options(scipen=999)  

meta= read.table("L1_L3_metainfo_2.0", header=F, sep = "\t")
table<- read.table("esxH_HAP", header=T)
colnames(table) <- c("n","G","H","reg")

table$H10=substr(table$H,10,10)

colnames(meta) =c("G","SAM","AC","year","c_i", "reg", "c_o", "L", "sub_L")

meta <- subset(meta, select = c("G","reg","c_i"))

t <- read.tree("RAxML_bestTree.L1_country")

rooted_t= root(t, outgroup = "G01705", resolve.root = TRUE)
rooted_t= drop.tip(rooted_t,"G01705")


tips=c()
rooted_t$edge.length=(rooted_t$edge.length/sum(rooted_t$edge.length))

res_SE_ai=c()
res_SE_m=c()
res_SA=c()
res_EA=c()
res_EA_noM=c()
res_SAf=c()
res_SAm=c()
res=c()



for (r in (1:10000)){
    tips=c()
    for (z in (1:8)){

        rnum=runif(1, 0, 1)
        tot=0
        for (i in (1:length(rooted_t$edge.length))){
    
            tot=tot+rooted_t$edge.length[i]
            if (rnum < tot){
                br=i
                break
            }
        }

        if (is.na(t$tip.label[br])){
            t1 <-extract.clade(t, br)
            tips = c(tips,t1$tip)
        }else{
            tips=c(tips,t$tip.label[br])
        } 
    }

    tt= unique(tips)

    sub <- subset(meta,meta$G %in% tt)
 
    

    count <- nrow(subset(sub,sub$reg =="S Africa"))
    res_SAf=c(res_SAf,count)

    count <- nrow(subset(sub,sub$reg =="S America"))
    res_SAm=c(res_SAm,count)
    
    count <- nrow(subset(sub,sub$reg =="SE Asia (islands)"))
    res_SE_ai=c(res_SE_ai,count)
    
    count <- nrow(subset(sub,sub$reg =="SE Asia (mainland)"))
    res_SE_m=c(res_SE_m,count)
    
    count <- nrow(subset(sub,sub$reg =="S and C Asia"))
    res_SA=c(res_SA,count)
   
    count <- nrow(subset(sub,sub$reg =="E Africa"))
    res_EA=c(res_EA,count)
 
    count=nrow(subset(sub,(sub$reg=="E Africa") & (sub$c_i!="Madagascar")))   
    res_EA_noM=c(res_EA_noM,count)
   
}


obs=nrow(subset(table,((table$H10=="V") & (table$reg=="SE Asia (islands)"))))
res_temp = c(res_SE_ai,obs)

pos <- match(obs,res_temp)


low<-(rank(res_temp,ties.method = c("first"))/length(res_temp))
high<-(rank(res_temp,ties.method = c("last"))/length(res_temp))
    
res=rbind(res,c("Southeast Asia (islands)",low[pos],high[pos],NA,obs))

obs=nrow(subset(table,((table$H10=="V") & (table$reg=="SE Asia (mainland)"))))
res_temp = c(res_SE_m,obs)

pos <- match(obs,res_temp)


low<-(rank(res_temp,ties.method = c("first"))/length(res_temp))
high<-(rank(res_temp,ties.method = c("last"))/length(res_temp))
    
res=rbind(res,c("Southeast Asia (mainland)",low[pos],high[pos],NA,obs))

obs=nrow(subset(table,((table$H10=="V") & (table$reg=="S and C Asia"))))
res_temp = c(res_SA,obs)

pos <- match(obs,res_temp)


low<-(rank(res_temp,ties.method = c("first"))/length(res_temp))
high<-(rank(res_temp,ties.method = c("last"))/length(res_temp))
    
res=rbind(res,c("South Asia",low[pos],high[pos],NA,obs))

obs=nrow(subset(table,((table$H10=="V") & (table$reg=="E Africa"))))
res_temp = c(res_EA,obs)

pos <- match(obs,res_temp)


low<-(rank(res_temp,ties.method = c("first"))/length(res_temp))
high<-(rank(res_temp,ties.method = c("last"))/length(res_temp))
    
res=rbind(res,c("East Africa",low[pos],high[pos],NA,obs))



obs=nrow(subset(table,((table$H10=="V") & (table$reg=="E Africa") & (table$c_i!="Madagascar"))))
res_temp = c(res_EA_noM,obs)

pos <- match(obs,res_temp)


low<-(rank(res_temp,ties.method = c("first"))/length(res_temp))
high<-(rank(res_temp,ties.method = c("last"))/length(res_temp))
    
res=rbind(res,c("East Africa (no Madagascar)",low[pos],high[pos],NA,obs))


obs=nrow(subset(table,((table$H10=="V") & (table$reg=="S Africa"))))
res_temp = c(res_SA,obs)

pos <- match(obs,res_temp)


low<-(rank(res_temp,ties.method = c("first"))/length(res_temp))
high<-(rank(res_temp,ties.method = c("last"))/length(res_temp))

res=rbind(res,c("South Africa",low[pos],high[pos],NA,obs))



obs=nrow(subset(table,((table$H10=="V") & (table$reg=="S America"))))
res_temp = c(res_SAm,obs)

pos <- match(obs,res_temp)


low<-(rank(res_temp,ties.method = c("first"))/length(res_temp))
high<-(rank(res_temp,ties.method = c("last"))/length(res_temp))

res=rbind(res,c("South America",low[pos],high[pos],NA,obs))

#res[,2] <- as.numeric(as.character(res[,2]))
#res[,3] <- as.numeric(as.character(res[,3]))




for (i in 1:nrow(res)){
    if (res[i,2] < 0.5) {res[i,4]<-res[i,3]}
    if (res[i,2] > 0.5) {res[i,4]<-res[i,2]}
}

colnames(res)<-c("area","low_ep","high_ep","empirical p-value", "count")

for (i in 1:nrow(res)){
  if (as.numeric(res[i,4]) > 0.5) {res[i,4]=(1-as.numeric(res[i,4]))}
}


res <- subset(res, select = c("area","empirical p-value"))

write.table(res,file="phyC_A10V_mut.txt",row.names=F)
