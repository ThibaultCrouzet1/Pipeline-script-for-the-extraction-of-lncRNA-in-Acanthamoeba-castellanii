install.packages("gplots")
install.packages("RColorBrewer")

library("gplots")
library("RColorBrewer")
library("heatmap.plus")
library("writexl")

##########################################################################################
#AC
##########################################################################################


Ac_lncRNA_expression_level <- read.delim(file = "~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/Ac/Ac_lncRNA_expression_level.tab"
                                    ,dec = ","  )#,row.names = 1)

Ac_df <- Ac_lncRNA_expression_level[,-24:-30]


df_module <- Ac_df #[Ac_df$ModuleID == "Module_1",]


df_module_mean <- data.frame(Transcript_ID = Ac_df$Transcript_ID , P0 = NA, P1 = NA , P4 = NA , P8 = NA, E1 = NA , E4 = NA , E8 = NA, ModuleID = Ac_df$ModuleID)# , row.names = 1 )
#to the means of the condition replicates
for (r in 1:length(row.names(df_module)) ){
  
  df_module_mean$P0[r] <- rowMeans(df_module[ r,2:4]  )
  df_module_mean$P1[r] <- rowMeans(df_module[ r,5:7] )
  df_module_mean$P4[r] <- rowMeans(df_module[ r,8:10] )
  df_module_mean$P8[r] <- rowMeans(df_module[ r,11:13] )
  df_module_mean$E1[r] <- rowMeans(df_module[ r,14:16] )
  df_module_mean$E4[r] <- rowMeans(df_module[ r,17:19] )
  df_module_mean$E8[r] <- rowMeans(df_module[ r,20:22] )
  df_module_mean$ModuleID[r] <- df_module$ModuleID[r]
  
}
write_xlsx(df_module_mean, "~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/Ac/Ac_lncRNA_expression_level_average_replicat.xlsx")


Ac_input <- data.matrix(df_module_mean[,1:7], rownames.force = NA)



Ac_condition_colors <- unlist(lapply(colnames(df_module_mean),function(x){
  if(grepl("P",x)) '#80FF00' #green
  else if(grepl('E',x)) '#FF8000' #orange
  
}))



Ac_groupe <- data.frame(module_name = as.factor(df_module_mean$ModuleID), module_nb = as.numeric(as.factor(df_module_mean$ModuleID)) )


rainbow_array = c("#0000FF","#006600" ,"#FF0033", "#FFFF00", "#9900FF", "#FF9933", "#000000")

Ac_module_color <- data.frame(module_name = Ac_groupe$module_name , color = "" )




for( i in 1:length(Ac_groupe$module_nb) )
{
  Ac_module_color$color[i] <- rainbow_array[Ac_groupe$module_nb[i]]
}





pdf(file = "~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/Ac/Ac_lncRNA_heatmap_replicate_average.pdf")

svg(file = "~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/Ac/Ac_lncRNA_heatmap_replicate_average.svg")

heatmap.2(Ac_input, trace="none", density="none", col=bluered(20), 
          cexRow=1, 
          cexCol=1, 
          margins = c(5,14),
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          scale="column",
          ColSideColors = Ac_condition_colors, 
          RowSideColors = Ac_module_color$color,
          Colv = NA,
          dendrogram = "row",
          key.xlab = "Expression Level"
          )


legend(0.8,1,legend = c("Condition P","Condition E"),fill=c('#80FF00','#FF8000'),cex=0.75)
legend(0,0.5,legend = unique(Ac_module_color$module_name),fill=unique(Ac_module_color$color),cex=0.75,bg = "white")

dev.off()





##########################################################################################
#ACPA
##########################################################################################

AcPa_lncRNA_expression_level <- read.delim("~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/AcPa/AcPa_lncRNA_expression_level.tab"
                                           ,dec = "," )

AcPa_df <- data.frame(AcPa_lncRNA_expression_level[,-18:-25],Taxonomy = AcPa_lncRNA_expression_level$Taxonomy)



df_module <- AcPa_df #[Ac_df$ModuleID == "Module_1",]


df_module_mean <- data.frame(Transcript_ID = df_module$TranscriptID , P0 = NA, P1 = NA , P4 = NA , E1 = NA , E4 = NA , ModuleID = df_module$ModuleID , Taxonomy = df_module$Taxonomy )#, row.names = 1 )
#to do the means of the condition replicates
for (r in 1:length(row.names(df_module)) ){
  
  df_module_mean$P0[r] <- rowMeans(df_module[r,8:10] )
  df_module_mean$P1[r] <- rowMeans(df_module[r,11:13] )
  df_module_mean$P4[r] <- rowMeans(df_module[r,14:16] )
  df_module_mean$E1[r] <- rowMeans(df_module[r,2:4]  )
  df_module_mean$E4[r] <- rowMeans(df_module[r,5:7] )
  df_module_mean$ModuleID[r] <- df_module$ModuleID[r]
  df_module_mean$Taxonomy[r] <- df_module$Taxonomy[r]
  
}
AcPa_input <- data.matrix(df_module_mean[,1:5], rownames.force = NA)


write_xlsx(df_module_mean, "~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/AcPa/AcPa_lncRNA_expression_level_average_replicat.xlsx")


AcPa_condition_colors <- unlist(lapply(colnames(df_module_mean[,1:5]),function(x){
  if(grepl("P",x)) '#80FF00' #green
  else if(grepl('E',x)) '#FF8000' #orange

}))



AcPa_groupe <- data.frame(module_name = as.factor(df_module_mean$ModuleID), module_nb = as.numeric(as.factor(df_module_mean$ModuleID)) )

AcPa_module_color <- data.frame(module_name = AcPa_groupe$module_name , color = "" )

AcPa_module_color$color <- brewer.pal(12,"Set3")[AcPa_groupe$module_nb]


AcPa_origine <- data.frame(origine_name = as.factor(df_module_mean$Taxonomy), origine_nb = as.numeric(as.factor(df_module_mean$Taxonomy)) )
AcPa_origine_color <- data.frame(origine_name = AcPa_origine$origine_name , color = "" )
AcPa_origine_color$color <- brewer.pal(3,"Set1")[AcPa_origine$origine_nb]
myCol <- cbind(AcPa_module_color$color, AcPa_origine_color$color)

pdf(file = "~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/AcPa/AcPa_lncRNA_heatmap.pdf")
svg(file = "~/boulo/cours/Cours-Gphy-M2/stage/fichier serveur/06.EXPRESSION_LEVEL/6.3.LNCRNA_EXPRESSION_LEVEL/AcPa/AcPa_lncRNA_replicat_average_heatmap.svg")

heatmap.plus(AcPa_input, 
             trace="none", 
             density="none", 
             col=bluered(20), 
          cexRow=1, 
          cexCol=1, 
          margins = c(25,35),
          hclust=function(x) hclust(x,method="complete"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          scale="column",
          ColSideColors = cbind(AcPa_condition_colors , NA ), 
          RowSideColors = cbind(AcPa_module_color$color, AcPa_origine_color$color) ,
          Colv = NA,
          Rowv = NA,
          main = "Expression Level"
          )


legend(0.8,0.65,legend = c("Condition P","Condition E"),fill=c('#80FF00','#FF8000'),cex=0.75,bg = "white")
legend(0,0.4,legend = c(unique(AcPa_module_color$module_name),unique(AcPa_origine_color$origine_name) ),fill=c(unique(AcPa_module_color$color),unique(AcPa_origine_color$color) ),cex=0.55,bg = "white")

dev.off()













