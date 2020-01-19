
###BioconductorTopGo

library(topGO)
library(ALL)
data(ALL)


# En esta etapa del analisis se crea un Objeto topGOdata
# que almacena la lista de genes con  valores de expresion diferencial para cada gen 
# y el mapa de anotaciones GO, que asocia cada anotacion con uno o mas genes.


## Anotaciones GO para Arabidopsis Thaliana

library(org.At.tair.db)

org.At.tairGO
columns(org.At.tairGO)
View(org.At.tair.db)

typeof(org.At.tair.db)
Go_terms_and_categoria = toTable(org.At.tairGO)
View(Go_terms_and_categoria)

# las anotaciones Go se dividen por ontologia; 
# MF (Molecular funciton), BP (Biological Proces) y CC (Celular component)

Go_terms_and_categoria_MF<- subset(Go_terms_and_categoria, Ontology=='MF')
View(Go_terms_and_categoria_BP)

Go_term_MF <- data.frame(Go_terms_and_categoria_MF$gene_id, Go_terms_and_categoria_MF$go_id)
View(Go_term_MF)

write.table(Go_term_MF, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/GO/Go_MF.txt", sep="\t")

# Crear un archivo con los GO para cada categoria ontologica. (Go_MF.txt, GO_BP.txt, GO_CC.txt) 



## LISTAS DE VECTORES

# Los genes y sus respectivos GO se presentan como una lista de vectores, 
#para eso, hay que convertir cada uno de los archivos con los GO por categoria en listas de vectores 

unilist2<-unique(Go_terms_and_categoria$gene_id)

GOlist<-vector("list",length(unilist2))

for (i in 1:length(unilist2)){
  idx<-which(unilist2[i]==Go_terms_and_categoria$gene_id)
  GOlist[[i]]<-Go_terms_and_categoria[idx,"go_id"]
  cat(i,sep="\n")
}

names(GOlist)<-unilist2
str(GOlist)

save
save(list= GOlist_BP, file="C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/GOlist/GOlist_BP.RData")
cat(capture.output(print(GOlist4), file="C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/GOlist/GOlist_BP.txt"))

#Repetir para GOlist_BP.txt, GOlist_MF.txt, GOlist_CC.txt


## LIST OF INTERESTING GENES

# Gen universe son todos los genes de microarreglo 
# (en este caso los que tienen sinonimo en la nueva version V6, y a su ves ortólogos en A. thaliana)


Universo = read.table("file:///C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Gen_universe/gene_universe_uniq_ortologos_arabidopsis.txt", header = T, sep="\t")

View(Universo)

as.character(Universo)

## MY INTERESTING GENES 

Gene_interes = read.table(
  "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/ortologos/down_regulated_pten/Ortologos_one_to_one_down_regulated_pten.txt", header = T, sep="\t")

View(Gene_interes)
myInterestedGenes2 <- sample(Gene_interes$Ortholog.gene.name, 125)
View(myInterestedGenes2)

geneList_downreg_pten <- factor(as.integer(Universo$Ortholog.gene.name%in%myInterestedGenes2 ))
names(geneList_downreg_pten) <- Universo$Ortholog.gene.name
str(geneList_downreg_pten)
print(geneList_downreg_pten)

cat(capture.output(print(geneList), file="C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Gen_universe/Universo_and_genes_interesantes/Int_gen_downreg_pten_one_to_one"))

# Mismo prosedimiento para Ortologos_one_to_one_up_regulated_pten.txt


### CONSTRUCCIÓN DEL OBJETO TOP GO 

library(topGO)

source("http://bioconductor.org/biocLite.R")
biocLite("topGO")

source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("GenomicFeatures", "AnnotationDbi"))


# Go object 

install.packages("installr")


library(installr)

updateR()


GOdata_BP_all_up_pten <- new("topGOdata", 
                             ontology = "BP",
                             allGenes = geneList_up_pten_all,
                             annot = annFUN.gene2GO, 
                             gene2GO = GOlist)
GOdata



selGenes <- sample(a, 10)
gs <- geneScore(GOdata, whichGenes = selGenes)
print(gs)

sg <- sigGenes(GOdata2)
str(sg)
numSigGenes(GOdata2)

graph(GOdata2)

ug <- usedGO(GOdata)
head(ug)

sel.terms <- sample(usedGO(GOdata), 10)
num.ann.genes <- countGenesInTerm(GOdata, sel.terms) ## the number of annotated genes
num.ann.genes
ann.genes <- genesInTerm(GOdata, sel.terms) ## get the annotations
head(ann.genes)


ann.score <- scoresInTerm(GOdata, sel.terms)
head(ann.score)
ann.score <- scoresInTerm(GOdata, sel.terms, use.names = TRUE)
head(ann.score)
termStat(GOdata, sel.terms)

