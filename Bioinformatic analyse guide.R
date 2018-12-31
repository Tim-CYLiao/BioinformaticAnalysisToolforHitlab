
#============================================================================================================================================================================#
#the first code to run to set up require package
remove(list=ls())

#all file better put under this working direction folder
setwd("Pathway")
source("Bioinformatic analyse.R")
#source("Bioinformatic analyse.R",encoding="utf-8")
instalmypackage()
#R will ask update package,chose some(type s and enter,then type y for each require)

#============================================================================================================================================================================#
#search yeast protein description on Cyclop database

#reportdirection is excel(.xlsx only) file direction
#file format require column name of protein list must be "ID"(ex:YDR054C),sheetname must have a word "hit".
#following are example code
reportdirection<-"2017-05-23-Analyse result.xlsx"
#function code:
a<-Sys.time()
Cyclopssearch(reportdirection)
b<-Sys.time()-a
#type b can show up process time,test result: 30 sec for 124 hits
b
#============================================================================================================================================================================#
#search  protein annotation (protein name,pathway,GO ID) on uniprot database.
#Saccharomyces cerevisiae organismid<-559292
#Ecoli K1 K12 organismid<-c(405955,83333)
organismid<-559292
#file format require column name of protein list must be "ID"(ex:YDR054C),sheetname must have a word "hit".
#chose input "ID" format:"Systematic.name"(ex:YDR054C) or "Gene.name"(ex:CDC34)
searchcolumn<-"Systematic.name"
#it require  "hitlab proteome library.xlsx" as protein name database
librarydirection<-"hitlab proteome library.xlsx"
reportdirection<-"2017-05-23-Analyse result.xlsx"
#adding ranking by intensity, excel file must have column called "Intensityranking",and rank number inside it,will not show up if FALSE
addintensityranking<-FALSE
#Chose protein involved GO ID source database 2 option for "Uniprot" (only result on uniprot)or "GO" (including uniprot, panther, SGD etc)
GOIDdb<-"GO"
#function code:
uniprottarget<-c("id","genes","protein%20names","comment(FUNCTION)","comment(PATHWAY)","interactor","comment(DISEASE)","comment(SUBCELLULAR%20LOCATION)","go(biological%20process)","go(cellular%20component)","go(molecular%20function)","go-id","organism-id")
uniprottargetcolnames<-c("Uniprot Entry","Gene names","Protein names","Function","Pathway","Interactor","Relate disease","Subcellular location","Uniprot BP","Uniprot CC","Uniprot MF","go-id")
a<-Sys.time()
Uniprotsearch(reportdirection,librarydirection,organismid,GOIDdb,searchcolumn,addintensityranking,uniprottarget,uniprottargetcolnames)
b<-Sys.time()-a
#type b can show up process time,test result: 13 sec for yeast 124 hits
b
#============================================================================================================================================================================#

#collect protein seq(".fasta") from uniprot, this offer file to search common motif on MEME GLAM2
#input report(reportdirection) must be result from Uniprotsearch function
reportdirection<-"found common motif protein-uniprot.xlsx"
#function code:


a<-Sys.time()
Uniprotseqsearch(reportdirection)
b<-Sys.time()-a
#type b can show up process time,test result: 1.5 min for 124 hits
b
#============================================================================================================================================================================#
#search GO data from bioconductor GO.db, it will show up GO information(GO ID,Term,Ontology and definection) and their protein hit


#input report(reportdirection) must be result from Uniprotsearch function
reportdirection<-"2017-05-23-Analyse result-uniprot.xlsx"

#can chose hit of protein name wanna to show up
#three kind protein name can chose,"Systematic.name" ,"Gene.name", "Uniprot.Entry"
showupname<-"Gene.name"
#function code:

a<-Sys.time()
GOsearch(reportdirection,showupname)
b<-Sys.time()-a
#type b can show up process time,test result: 1 sec for 124 hits
b
#============================================================================================================================================================================#

#input report(reportdirection) must be result from GOsearch function
GOreportdir<-"2017-05-23-Analyse result-GO.xlsx"
#ancestorlayer is a factor for number of ancestor for each class,bigger number will have detail-class but also let classification more complicate 
ancestorlayercount<-9
#if searchedinputtarget is specific vector of GO ID ,it will search which children GOID(from protein) and protein is belong to it, and will close ancestorlayer function
#also,must give searchedinputtarget correct searchedinputtargetontologytype(ex:"BP","CC","MF").Remember, only can search 1 kind of ontology 1 time.
searchedinputtarget<-""
#searchedinputtargetontologytype is mean classify protein depend on previous GOID search result
searchedinputtargetontologytype<-"auto"

#function code:

a<-Sys.time()
for (x in 1:ancestorlayercount)
   {
   	   ancestorlayer<-x

   	   GOclassify(GOreportdir,uniprotreportdir,ancestorlayer,searchedinputtarget,searchedinputtargetontologytype)
   }

b<-Sys.time()-a
#type b can show up process time,test result: 8 sec for 124 hits
b
#============================================================================================================================================================================#

aaa

for(x in 2:nrow(aaa))
   {
   	   b<-grep(as.character(aaa[x,1]),as.character(organismidDBlist[[organismidDBcount]][,2]),value=T)
   	   if(length(b)==0)
   	      {
   	      	aaa[x,2]<-"Nss"

   	      }else
   	         {
   	         	if(length(b)==1)
   	         	   {
   	         	   	aaa[x,2]<-b
   	         	   	}else
   	         	   	   {
   	         	   	   	aaa[x,2]<-"M"
   	         	   	   }
   	         }
   }




   overlapGO<-read.csv("2017-05-23-Analyse result-GO classification with 3 ancestor layer.csv")
overlapGO[,"Protein.list"]
listresultlist<-vector("list",length=nrow(overlapGO))

   	   listresultlist<-strsplit(overlapGO[,"Protein.list"],";")
names(listresultlist)<-overlapGO[,1]
intersect()
overlapcc<-intersect(listresultlist[[8]],listresultlist[[9]])

overlapbpmb<-intersect (intersect(intersect(listresultlist[[3]],listresultlist[[4]]),listresultlist[[5]]),listresultlist[[7]])
overlapbp<-intersect(listresultlist[[1]],listresultlist[[2]])
for(x in 1:7)
   {
   	   overlapbp<-intersect(listresultlist[[x]],overlapbp)
   }
totalcc<-union(listresultlist[[8]],listresultlist[[9]])

totalmb<-union (union(union(listresultlist[[3]],listresultlist[[4]]),listresultlist[[5]]),listresultlist[[7]])
totalbp<-union(listresultlist[[1]],listresultlist[[2]])
for(x in 1:7)
   {
   	   totalbp<-union(listresultlist[[x]],totalbp)
   }







   overlapbpother<-intersect(intersect(listresultlist[[1]],listresultlist[[2]]),listresultlist[[6]])


   overlapcellularprocessi<-intersect(intersect(listresultlist[[1]],listresultlist[[4]]),listresultlist[[6]])
   overlapcellularprocessu<-union(union(listresultlist[[1]],listresultlist[[4]]),listresultlist[[6]])







