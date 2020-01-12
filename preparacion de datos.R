###BIOCONDUCTOR 

source("http://www.bioconductor.org/biocLite.R")

biocLite()
biocValid()

#https://www.bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf

source("http://bioconductor.org/biocLite.R")
biocLite("topGO")


# El punto de partida es una lista de genes y sus p-valores de expresion diferencial (si es la intencion ver expresion diferencial) 

## T-TEST: Analisis estadistico para valuar el cambio de expresion(fold_change), por gen, entre la linea wild type(wt_a, wt_b, wt_c) 
## y la linea mutada para el gen pten (ko_a, ko_b, ko_c)
# Se obtiene un p-valor por gen genes que indica si la diferencia de expresion entre las dos condiciones, por gen, es significativa (t test)

#https://www.bioconductor.org/help/course-materials/2015/Uruguay2015/day5-data_analysis.html

getwd()
setwd("C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/")

pten = read.table("Analysis_pten_entero.txt", sep = "\t", header = T, dec = ",")

#transformar columna de genes en nombre de filas 
dimnames(pten)[[1]] <- pten[,1]
pten = pten[,-1]

dim(pten)
View(pten)
colnames(pten)


# t test entre las dos condiciones (wt y ko) del microarreglo 

y <- t.test(pten[1,1:3], pten[1, 4:6])

#Función "t-test" para calcular la diferencia por gen entre el wt(control) y ko

ttestPten <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}

rawpvalue = apply(pten, 1, ttestPten, grp1 = c(1:3), grp2 = c(4:6))   
View (rawpvalue)

hist(rawpvalue)

#Log2 de los datos, calcular la media de cada gen por grupo. 
#Luego calcular el fold change entre los grupos (control versus dieta cetogénica). 


##log2.
pten1 = log2(pten)

#Media del control (wt)
control = apply(pten1[,1:3], 1, mean)

#Media del tratamiento (ko)
test = apply(pten1[, 4:6], 1, mean) 

#confirmar que tenemos un vector de numeros
class(control) 
class(test)

#Nuestros datos ya están transformados en log2; podemos calcular la diferencia entre las medias. 
#Este es nuestro log2 Fold Change o log2 Ratio == log2 (control / test)

foldchange <- control - test 

hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

#agregar p valor como columna en data frame rat 

pten1$pvalue <- rawpvalue
head(pten1)

pten1$foldchange <- foldchange
head(pten1)

pten1$fold_change <- NULL
head(pten1)

#Transform the p-value (-1*log(p-value)) and create a volcano plot using ggplot2.

results = cbind(foldchange, rawpvalue)
results = as.data.frame(results)
results$probename <- rownames(results)

library(ggplot2)
volcano = ggplot(data = results, aes(x = foldchange, y = -1*log10(rawpvalue)))
volcano + geom_point()


write.table(pten1, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/all_genes_and_pvalues_pten.txt", sep="\t")
