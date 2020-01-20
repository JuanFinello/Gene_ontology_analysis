
###BioconductorTopGo

if(!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install(c("topGO","ALL"))

library(topGO)
library(ALL)
data(ALL)


# En esta etapa del analisis se crea un Objeto topGOdata
# que almacena la lista de genes y el mapa de anotaciones GO, 
# que asocia cada anotacion con uno o mas genes.


## Anotaciones GO para Arabidopsis Thaliana

source("https://bioconductor.org/biocLite.R")
biocLite("bit")

library(org.At.tair.db)

org.At.tairGO
columns(org.At.tairGO)
View(org.At.tair.db)

typeof(org.At.tair.db)
Go_terms_and_categoria = toTable(org.At.tairGO)
View(Go_terms_and_categoria)

# las anotaciones Go se dividen por ontologia; 
# MF (Molecular funciton), BP (Biological Proces) y CC (Celular component)

Go_terms_and_categoria_BP<- subset(Go_terms_and_categoria, Ontology=='BP')
Go_term_BP <- data.frame(Go_terms_and_categoria_BP$gene_id, Go_terms_and_categoria_BP$go_id)
colnames(Go_term_BP)<- c("gene_id", "go_id")


write.table(Go_term_BP, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/GO/Go_BP.txt", sep="\t")

# Crear un archivo con los GO para cada categoria ontologica. (Go_MF.txt, GO_BP.txt, GO_CC.txt) 



## LISTAS DE VECTORES

# Los genes y sus respectivos GO se presentan como una lista de vectores, 
#para eso, hay que convertir cada uno de los archivos con los GO por categoria ontologica en listas de vectores 

unilist2<-unique(Go_term_BP$gene_id)

GOlist_BP<-vector("list",length(unilist2))

for (i in 1:length(unilist2)){
  idx<-which(unilist2[i]==Go_terms_and_categoria_BP$gene_id)
  GOlist_BP[[i]]<-Go_terms_and_categoria_BP[idx,"go_id"]
  cat(i,sep="\n")
}

names(GOlist_BP)<-unilist2
str(GOlist_BP)

#save

cat(capture.output(print(GOlist_BP), file="C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/GOlist/GOlist_BP.txt"))

#Repetir para GOlist_BP.txt, GOlist_MF.txt, GOlist_CC.txt


## LIST OF INTERESTING GENES

# Gen universe son todos los genes de microarreglo 
# (en este caso los que tienen sinonimo, en la nueva version V6, de anotaciones de genes de P.patens,
# y a su ves ort?logos en A. thaliana)


Universo = read.table("file:///C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Gen_universe/gene_universe_uniq_ortologos_arabidopsis.txt", header = T, sep="\t")

View(Universo)

as.character(Universo)

## MY INTERESTING GENES 

# El paso de versiones de anotación (V3 a V6) y luego a ortologos hace que perdamos los p-valores para cada secuencia, 
# En este escenario, cuando solo se proporciona una lista de genes interesantes (en su condicion up or down), 
# solo puede usarse estadísticas de pruebas basadas en recuentos de genes, como el de Fisher.

Gene_interes_down = read.table(
  "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/ortologos/down_regulated_pten/Ortologos_one_to_one_down_regulated_pten.txt", header = T, sep="\t")

View(Gene_interes_down)
myInterestedGenes2 <- sample(Gene_interes_down$Ortholog.gene.name, 249)
View(myInterestedGenes2)


#The geneList object is a named factor that indicates which genes are interesting and which not

geneList_downreg_pten <- factor(as.integer(Universo$Ortholog.gene.name%in%myInterestedGenes2 ))
names(geneList_downreg_pten) <- Universo$Ortholog.gene.name
str(geneList_downreg_pten)
print(geneList_downreg_pten)

cat(capture.output(print(geneList_downreg_pten), file="C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Gen_universe/Universo_and_genes_interesantes/Int_gen_downreg_pten_one_to_one"))

# Mismo procedimiento para Ortologos Up regulated

Gene_interes_up = read.table(
  "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/ortologos/up_regulated_pten/Ortologos_all_up_reg_pten.txt", header = T, sep="\t")

View(Gene_interes_up)
myInterestedGenes1 <- sample(Gene_interes_up$Ortholog.gene.name, 353)
View(myInterestedGenes1)

geneList_upreg_pten <- factor(as.integer(Universo$Ortholog.gene.name%in%myInterestedGenes1 ))
names(geneList_upreg_pten) <- Universo$Ortholog.gene.name
str(geneList_upreg_pten)
print(geneList_upreg_pten)

cat(capture.output(print(geneList_upreg_pten), file="C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Gen_universe/Universo_and_genes_interesantes/Int_gen_upreg_pten_one_to_one"))




### CONSTRUCCI?N DEL OBJETO TOP GO 

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

# Go_Up

# Para construir el onjeto topGOdata la lista d egenes interesante, la ontologia de interes, y las anotaciones GO  


GOdata_BP_all_up_pten <- new("topGOdata", 
                             ontology = "BP",
                             allGenes = geneList_up_pten_all,
                             annot = annFUN.gene2GO, 
                             gene2GO = GOlist_BP)

getClass(Class, where = topenv(parent.frame()))

stop(gettextf("%s is not a defined class", dQuote(Class)), domain = NA)



