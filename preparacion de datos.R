### PREPARACIÓN DE DATOS PARA ANALISIS GO

# Para el analisis GO, El punto de partida es una lista de genes y sus p-valores de 
# expresion diferencial (si es la intencion ver expresion diferencial) 
# El resultado del microarray es una trabla con 32559 genes, 
# y sus valores de expresion en dos condiciones, wt y ko, control y tratamiento respectivamente 
# con tres repeticiones para cada caso.



## T-TEST:

# Analisis estadistico para valuar el cambio de expresion (fold_change), por gen, entre la linea wild type (wt) 
# y la linea mutada para el gen pten (ko)
# Se obtiene un p-valor por gen genes que indica si la diferencia de expresion entre las dos condiciones, por gen, es significativa (t test)

# fuente https://www.bioconductor.org/help/course-materials/2015/Uruguay2015/day5-data_analysis.html

getwd()
setwd("C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/")

# Importar archivo con resultados de microarray

pten = read.table("Analysis_pten_entero.txt", sep = "\t", header = T, dec = ",")

#transformar columna de genes en nombre de filas 

dimnames(pten)[[1]] <- pten[,1]
pten = pten[,-1]

dim(pten)
View(pten)
colnames(pten)


# t test entre las dos condiciones (wt y ko) del microarreglo 

y <- t.test(pten[1,1:3], pten[1, 4:6])

# Función "t-test" 

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

# Log2 de los datos, calcular la media de cada gen por grupo. 
# Luego calcular el fold change entre los grupos (control versus dieta cetogénica). 



## log2.

# Log2 de los datos, calcular la media de cada gen por grupo. 
# Luego calcular el fold change entre los grupos (control versus dieta cetogénica). 


pten1 = log2(pten)

# Media del control (wt)
control = apply(pten1[,1:3], 1, mean)

# Media del tratamiento (ko)
test = apply(pten1[, 4:6], 1, mean) 

# confirmar que tenemos un vector de numeros
class(control) 
class(test)

# Nuestros datos ya están transformados en log2; podemos calcular la diferencia entre las medias. 
# Este es nuestro log2 Fold Change o log2 Ratio == log2 (control / test)

foldchange <- control - test 

hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

# Agregar p valor como columna en data frame rat 

pten1$pvalue <- rawpvalue
head(pten1)

pten1$foldchange <- foldchange
head(pten1)

pten1$fold_change <- NULL
head(pten1)

# Transform the p-value (-1*log(p-value)) and create a volcano plot using ggplot2.

results = cbind(foldchange, rawpvalue)
results = as.data.frame(results)
results$probename <- rownames(results)

library(ggplot2)
volcano = ggplot(data = results, aes(x = foldchange, y = -1*log10(rawpvalue)))
volcano + geom_point()


write.table(pten1, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/all_genes_and_pvalues_pten.txt", sep="\t")


### Filtrar por p valor significativo (Menor a 0.05)  

source(	broom_0.5.0.tar.gzm)

install.packages("broom")

library(broom)

# Now add row names as colum

library(tibble)
pten1<- rownames_to_column(pten1, var = "Genes")
head(pten1)
View(pten1)


##Filtrar genes segun tengan p-valor significativo; menor que 0.05

seleccion_pten_pvalor1<- as.data.frame(pten1$Genes[pten1$pvalue < 0.05])
seleccion_pten_pvalor2<- as.data.frame(pten1$pvalue[pten1$pvalue < 0.05])
seleccion_pten_pvalor3<- as.data.frame(pten1$foldchange[pten1$pvalue < 0.05])

head(seleccion_pten_pvalor3)
View(seleccion_pten_pvalor3)

class(seleccion_pten_pvalor1)

#crear nuevo Data.frame con valores filtrados 

filtrado_por_pvalue<- data.frame(seleccion_pten_pvalor1,seleccion_pten_pvalor2, seleccion_pten_pvalor3)
head(filtrado_por_pvalue)
View(filtrado_por_pvalue)

colnames(filtrado_por_pvalue) <- c("Genes", "pvalue", "foldchange")

write.table(filtrado_por_pvalue, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Pten_filtrado_por_pval_signif.txt", sep="\t")



## Clasificar por fold change, segun sea positivo (Up regulated) o negativo (Down regulated)
# Explicacion del valor de corte del fold change (+0.5849 para Up y - 0.5849)  https://www.biostars.org/p/101727/

#positivo "Up_regulated"

seleccion_pten_foldchange1<- as.data.frame(filtrado_por_pvalue$Genes[filtrado_por_pvalue$foldchange> 0.5849])
seleccion_pten_foldchange2<- as.data.frame(filtrado_por_pvalue$pvalue[filtrado_por_pvalue$foldchange > 0.5849])
seleccion_pten_foldchange3<- as.data.frame(filtrado_por_pvalue$foldchange[filtrado_por_pvalue$foldchange > 0.58498])


head(seleccion_pten_foldchange3)
View(seleccion_pten_foldchange2)

Up_regulated<- data.frame(seleccion_pten_foldchange1,seleccion_pten_foldchange2, seleccion_pten_foldchange3)
colnames(Up_regulated) <- c("Genes", "pvalue", "foldchange")

head(Up_regulated)
View(Up_regulated)


#crear un archivo de texto separado por tab con el data frame de los genes upregulated

write.table(Up_regulated, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Up_regulated_pten.txt", sep="\t")





## Negativo "Down_regulated"

seleccion_pten_foldchange4<- as.data.frame(filtrado_por_pvalue$Genes[filtrado_por_pvalue$foldchange < -0.5849])
seleccion_pten_foldchange5<- as.data.frame(filtrado_por_pvalue$pvalue[filtrado_por_pvalue$foldchange < -0.5849])
seleccion_pten_foldchange6<- as.data.frame(filtrado_por_pvalue$foldchange[filtrado_por_pvalue$foldchange < -0.58498])


head(seleccion_pten_foldchange3)
View(seleccion_pten_foldchange2)

Down_regulated<- data.frame(seleccion_pten_foldchange4,seleccion_pten_foldchange5, seleccion_pten_foldchange6)
colnames(Down_regulated) <- c("Genes", "pvalue", "foldchange")
head(Down_regulated)
View(Down_regulated)


#crear un archivo de texto de los genes downregulated

write.table(Down_regulated, "C:/Users/Admin/Desktop/tesina/bioinformatica microarreglo/archivos definitivos/Down_regulated_pten.txt", sep="\t")



