
### Analisis de enriquecimiento 
## Correr el analisis de enriquecimiento con los objetos GO construidos en el paso anterior 


library(topGO)
resultFis_CC_up_pten_all <- runTest(Godata_uf, algorithm = "classic", statistic = "fisher")
??runTest

##Interpretation and visualization of results

head(score(resultFis_CC_up_pten_all))

pvalFis <- score(resultFis_CC_up_pten_all)
head(pvalFis)
hist(pvalFis, 50, xlab = "p-values")


##Summarising the results

allRes <- GenTable(GOdata_CC_up_pten_all, classic = resultFis_CC_up_pten_all, ranksOf = "classic", topNodes = 20)
allRes


###Visualising the GO structure

source("https://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")

library(Rgraphviz)

??Rgraphviz



#Resultados en formato tabla
weight_fisher_result=runTest(GOdata_BP_up_pten_all, algorithm='weight01', statistic='fisher')
allGO=usedGO(GOdata_BP_up_pten_all)
allRess <- GenTable(GOdata_BP_up_pten_all, weightFisher= weight_fisher_result, orderBy = "weightFisher", ranksOf = "classicFisher", topNodes=20)

View(allRess)

write.table(allRes, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/resultados/up_all_pten/BP/tabla_up_pten_all_BP.txt", sep="\t")


par(cex = 0.2)

showSigOfNodes(GOdata_BP_up_pten_all, score(weight_fisher_result), firstSigNodes = 5, useInfo = 'all', useFullNames = TRUE, .NO.CHAR=100)

printGraph(GOdata_BP_up_pten_all, resultFis_BP_up_pten_all, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)



#performing BH correction on our p values
p.adj=round(p.adjust(allRess$weightFisher,method="BH"),digits = 4)

# create the file with all the statistics from GO analysis
all_res_final=cbind(allRess,p.adj)
all_res_final=all_res_final[order(all_res_final$p.adj),]
View(all_res_final)


#get list of significant GO after multiple testing correction
results.table.bh=all_res_final[which(all_res_final$p.adj<=0.05),]
View(results.table.bh)

#### para buscar por go especifico 

allGO["GO:0007623"]
RetrivedGenes <- lapply(allGO,function(x) x[x %in% myInterestedGenes3] )
RetrivedGenes[["GO:0007623"]]




###Get all the genes in your significant GO TERMS
# https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/

myterms = results.table.bh$GO.ID # change it to results.table.bh$GO.ID if working with BH corrected values
mygenes = genesInTerm(GOdata_BP_up_pten_all, myterms)

var=c()
for (i in 1:length(myterms))
{
  myterm=myterms[i]
  mygenesforterm= mygenes[myterm][[1]]
  mygenesforterm=paste(mygenesforterm, collapse=',')
  var[i]=paste("GOTerm",myterm,"genes-",mygenesforterm)
}
View(var)
write.table(var,"C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/genetoGOmapping.txt",sep="\t",quote=F)




