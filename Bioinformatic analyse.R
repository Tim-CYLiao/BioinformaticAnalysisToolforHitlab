#============================================================================================================================================================================#
instalmypackage<-function()
   {
       
       Pkg <- c("fields","spam","maps","base","openxlsx","xts","stringr","jpeg","RCurl","XML","dplyr","tseries")
       inst <- Pkg %in% installed.packages()
       if(length(Pkg[!inst]) > 0) install.packages(Pkg[!inst])
       instpackages <- lapply(Pkg, library, character.only=TRUE)
       source("https://bioconductor.org/biocLite.R")
       biocLite("GO.db")
       library("GO.db")
       #library("GO.db")


   }
#============================================================================================================================================================================#
username<-function()
   {
       sysinformation<-Sys.getenv()
          if(length(grep("Apple",names(Sys.getenv(names=TRUE))))==1)
             {
                 systemname<-sysinformation[["USER"]]
             }else
                {
                    systemname<-sysinformation[["USERNAME"]]
                }
       systemname<-as.character(systemname)
       return(systemname)
   } 
#============================================================================================================================================================================#

informationlistcombine<-function(inputlist,splitunit)
   {
   	  outputvector<-vector("character",length(inputlist))
   	  for(x in 1:length(inputlist))
   	     {
   	     	if(length(inputlist[[x]])!=0)
   	     	   {
   	     	   	   for(listunitvectercount in 1:length(inputlist[[x]]))
   	     	   	      {
   	     	   	      	   if(listunitvectercount==1)
   	     	   	      	      {
   	     	   	      	      	tempunit<-inputlist[[x]][listunitvectercount]
   	     	   	      	      }else
   	     	   	      	         {
   	     	   	      	         	tempunit<-paste(tempunit,inputlist[[x]][listunitvectercount],sep=splitunit)
   	     	   	      	         }
   	     	   	      }
   	     	   	   outputvector[x]<-tempunit
   	     	   }
   	     	
   	     }
   	  return(outputvector)
   } 
#============================================================================================================================================================================#
informationcombine<-function(input,splitunit)
   {
      outputvector<-""
      if(length(input)!=0)
             {
                 for(listunitvectercount in 1:length(input))
                    {
                         if(listunitvectercount==1)
                            {
                              tempunit<-input[listunitvectercount]
                            }else
                               {
                                tempunit<-paste(tempunit,input[listunitvectercount],sep=splitunit)
                               }
                    }
                 outputvector<-tempunit
             }

      return(outputvector)
   } 
#============================================================================================================================================================================#

Cyclopssearch<-function(reportdirection)
   {
       #Download genelist first
       CycloDB<-"http://cyclops.ccbr.utoronto.ca/genelist.html"
       Cyclopshtml <- getURL(CycloDB,cainfo=system.file("CurlSSL","cacert.pem",package="RCurl"))
       temp<-strsplit(Cyclopshtml,"<tr align=center bgcolor=#CCCCCC><td><font class=words><a href=")[[1]]
       temp<-temp[2:length(temp)]
       DBlist<-matrix(nrow=length(temp),ncol=4)
       title<-c("ID","Common name","Alias","Description")
       colnames(DBlist)<-title
       rownames(DBlist)<-c(seq(1:length(temp)))
       for(x in 1:length(temp))
          {
               temp02<-strsplit(temp[x],"<font class=words>")[[1]]
               temp02<-gsub("</font></td><td>",replacement="",temp02)
               temp02[4]<-strsplit(temp02[4],"</font></td></tr>\n")[[1]][1]
               temp03<-strsplit(temp02[1],"<")[[1]][1]
               temp02[1]<-strsplit(temp03,">")[[1]][2]
               DBlist[x,]<-temp02

          }

       report<-loadWorkbook(reportdirection)
       assayanalysisrsult<-createWorkbook(creator =username())
       reportnames<-grep("hit",names(report),value=T)
       analyseIDresult<-vector("list",length(names(report)))
       for(filecount in 1:length(reportnames))
          {
             #讀取檔案，與製作新的表格
             analyseIDresult[[filecount]]<-readWorkbook(report, sheet = reportnames[filecount],startRow = 1,colNames = TRUE,rowNames = FALSE)
             analyseIDresult[[filecount]]<-merge(analyseIDresult[[filecount]],DBlist,all.x=T,by="ID")
             missplace<-which(is.na(analyseIDresult[[filecount]]$Description))
             for(colcount in 1:length(colnames(analyseIDresult[[filecount]])))
                {
                    analyseIDresult[[filecount]][,colcount]<-as.character(analyseIDresult[[filecount]][,colcount])
                }
             if(length(missplace)!=0)
                {
                    for(placecount in 1:length(missplace))
                       {
                           if(grepl("*",analyseIDresult[[filecount]][missplace[placecount],"ID"])==TRUE)
                              {
                                   proteinID<-gsub("*",replacement="",analyseIDresult[[filecount]][missplace[placecount],"ID"])
                              }else
                                 {
                                     proteinID<-analyseIDresult[[filecount]][missplace[placecount],"ID"]
                                 }
                           
                           cyclopspath<-paste("http://cyclops.ccbr.utoronto.ca/SEARCH/search.pl?gene=",proteinID,sep="")
                           Cyclopshtml <- getURL(cyclopspath,cainfo=system.file("CurlSSL","cacert.pem",package="RCurl"))
                           temp<-strsplit(Cyclopshtml,"class=words>")[[1]][3:6]
                           IDsearchresult<-vector(length=length(temp))
                           for(x in 1:length(temp))
                              {
                                   IDsearchresult[x]<-strsplit(strsplit(temp[x],"<b>")[[1]][2],"</b>")[[1]][1]
                              }
                           IDsearchresult[1]<-proteinID
                           IDsearchresult<-as.matrix(t(IDsearchresult))
                           colnames(IDsearchresult)<-title
                           for(y in 2:length(title))
                              {
                                   
                                   analyseIDresult[[filecount]][missplace[placecount],title[y]]<-IDsearchresult[1,title[y]]
                              }


                           rm(IDsearchresult)
                       
                       }
                }

              addWorksheet(assayanalysisrsult, reportnames[filecount])
              writeDataTable(assayanalysisrsult, reportnames[filecount], analyseIDresult[[filecount]], startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = FALSE,tableStyle = "TableStyleLight9",withFilter = FALSE)

          }

          reportdirection2<-gsub("-Cyclop",replacement="",reportdirection)
          reportdirection2<-gsub(".xlsx",replacement="",reportdirection2)
          reportdirection2<-paste(reportdirection,"-","Cyclop",".xlsx",sep="")
          saveWorkbook(assayanalysisrsult, reportdirection2, overwrite = TRUE)
   }

#============================================================================================================================================================================#
#Saccharomyces cerevisiae organismid<-559292
#Ecoli K1 K12 organismid<-c(405955,83333)
#uniprotsearch(reportdirection,organismid)
searchcolumn<-"Systematic.name"
#librarydirection<-
#reportdirection<-
addintensityranking<-TRUE
#Chose protein involved GO ID source database 2 option for "Uniprot" (only result on uniprot)or "GO" (including uniprot, panther, SGD etc)


GOIDdb<-"GO"
#uniprottarget<-c("id","genes","protein%20names","pathway","go-id","organism-id")
#uniprottargetcolnames<-c("Uniprot Entry","Gene names","Protein names","Pathway","go-id")
Uniprotsearch<-function(reportdirection,librarydirection,organismid,GOIDdb,searchcolumn,addintensityranking,uniprottarget,uniprottargetcolnames)
   {

       organismidDBlist<-vector("list",length=length(organismid))
       for(organismidDBcount in 1:length(organismid))
          {
               #uniprottarget<-c("id","genes","protein%20names","pathway","go","organism-id")
               uniprottargetbind<-uniprottarget[1]
               for(bindcount in 2:length(uniprottarget))
                  {
                       uniprottargetbind<-paste(uniprottargetbind,uniprottarget[bindcount],sep=",")
                  }
               
               uniprotpath<-paste("http://www.uniprot.org/uniprot/?query=",organismid[organismidDBcount],"&format=tab&columns=",uniprottargetbind,sep="")
               uniprothtml <- getURL(uniprotpath,cainfo=system.file("CurlSSL","cacert.pem",package="RCurl"))
               temp<-strsplit(uniprothtml,"\n")[[1]]
               temp2<-grep(organismid[organismidDBcount],temp,value=T)
               searchlist<-matrix("Notfound",ncol=(length(uniprottarget)-1),nrow=length(temp2))
               searchlist<-as.data.frame(searchlist)
               #colnames(searchlist)<-c("Uniprot Entry","Gene names","Protein names","Pathway","Gene ontology (GO)")
               colnames(searchlist)<-uniprottargetcolnames
               for(colcount in 1:length(colnames(searchlist)))
                  {
                      searchlist[,colcount]<-as.character(searchlist[,colcount])
                  }
               
               for(dbproteincount in 1:length(temp2))
                  {
                      searchlist[dbproteincount,1:(length(uniprottarget)-1)]<-t(strsplit(temp2[dbproteincount],"\t")[[1]])[,1:(length(uniprottarget)-1)]
                  }

               if(length(intersect(uniprottargetcolnames,"Function"))==1)
                  {
                      tempfunction<-gsub("FUNCTION: ",replacement="",searchlist[,"Function"])
                      tempfunctionlist<-strsplit(tempfunction,"; ")
                      for(tempfunctionlistcount in 1:length(tempfunctionlist))
                         {
                             tempfunctionlist[[tempfunctionlistcount]]<-sapply(strsplit(tempfunctionlist[[tempfunctionlistcount]],as.character("ECO")),"[[",1)
                                     
                         }
                      splitunit<-";"
                      searchlist[,"Function"]<-informationlistcombine(tempfunctionlist,splitunit)

                  }
               if(length(intersect(uniprottargetcolnames,"Subcellular location"))==1)
                  {
                      tempfunction<-gsub("SUBCELLULAR LOCATION: ",replacement="",searchlist[,"Subcellular location"])
                      tempfunctionlist<-strsplit(tempfunction,"; ")
                      for(tempfunctionlistcount in 1:length(tempfunctionlist))
                         {
                             tempfunctionlist[[tempfunctionlistcount]]<-sapply(strsplit(tempfunctionlist[[tempfunctionlistcount]],as.character("ECO")),"[[",1)
                                     
                         }
                      splitunit<-";"
                      searchlist[,"Subcellular location"]<-informationlistcombine(tempfunctionlist,splitunit)
                  }
               if(length(intersect(uniprottargetcolnames,"Pathway"))==1)
                  {
                      searchlist[,"Pathway"]<-gsub("PATHWAY: ",replacement="",searchlist[,"Pathway"])
                  }
               organismidDBlist[[organismidDBcount]]<-searchlist
          }
       report<-loadWorkbook(reportdirection)
       libraryname<-loadWorkbook(librarydirection)
       librarydatalist<-vector("list",(length(organismid)+1))
       #合併兩個library
       for(organismidDBcount in 1:length(organismid))
          {
              librarysheetname<-grep(organismid[organismidDBcount],names(libraryname),value=T)
              librarydatalist[[organismidDBcount]]<-readWorkbook(libraryname, sheet = librarysheetname,startRow = 1,colNames = TRUE,rowNames = FALSE)
          	  if(organismidDBcount==1)
          	     {
          	     	  librarydata<-librarydatalist[[organismidDBcount]]
          	     }else
          	        {
          	        	 librarydata<-rbind(librarydata,librarydatalist[[organismidDBcount]])
          	        }
          }
       

       assayanalysisrsult<-createWorkbook(creator =username())
       
       reportnames<-grep("hit",names(report),value=T)
       analyseIDresult<-vector("list",length(names(report)))
       for(filecount in 1:length(reportnames))
          {
             #讀取檔案，與製作新的表格
             addWorksheet(assayanalysisrsult, reportnames[filecount])
             analyseIDresult[[filecount]]<-readWorkbook(report, sheet = reportnames[filecount],startRow = 1,colNames = TRUE,rowNames = FALSE)
             for(colcount in 1:length(colnames(analyseIDresult[[filecount]])))
                {
                    analyseIDresult[[filecount]][,colcount]<-as.character(analyseIDresult[[filecount]][,colcount])
                }
             searchlist<-matrix("Notfound",ncol=1,nrow=nrow(analyseIDresult[[filecount]]))
             #searchlist[,4:6]<-""
             searchlist<-as.data.frame(searchlist)
             for(colcount in 1:length(colnames(searchlist)))
                {
                    searchlist[,colcount]<-as.character(searchlist[,colcount])
                }
             #colnames(searchlist)<-c("ID","JWID","Gene names","Protein names","Pathway","Gene ontology (GO)")
             #colnames(searchlist)<-c("Uniprot Entry","Gene names","Protein names","Pathway","go-id")
             
             colnames(searchlist)<-c(searchcolumn)
             searchlist[,searchcolumn]<-analyseIDresult[[filecount]][,"ID"]
  
             if(searchcolumn=="Systematic.name")
                {
                	 searcheddatacolarrange<-c(1,2,4)
                	 searchcolumn2<-"Systematic name"
                }else
                   {
                   	   searcheddatacolarrange<-c(2,1,4)
                   	   searchcolumn2<-"Gene name"
                   } 
         
             searcheddata<-merge(searchlist,librarydatalist[[organismidDBcount]],all.x=TRUE,by=searchcolumn)[,searcheddatacolarrange]
             
             colnames(searcheddata)<-c("Systematic name","Gene name","Uniprot Entry")

             outputresult<-merge(searcheddata,organismidDBlist[[organismidDBcount]],all.x=TRUE,by="Uniprot Entry")
             outputresult<-outputresult[,c(2,3,5,1,c(6:ncol(outputresult)))]
             

             if(GOIDdb=="GO")
                {
                    templookingentryid<-grep("Not found",outputresult[,"Uniprot Entry"],value=TRUE,invert=TRUE)
                    templookingentryid<-grep("More than 1",templookingentryid,value=TRUE,invert=TRUE)
                    uniprotpath<-paste("http://www.ebi.ac.uk/QuickGO/GAnnotation?protein=",informationcombine(templookingentryid,","),"&format=tsv&col=proteinID,goID&limit=-1",sep="")
                    uniprothtml <- getURL(uniprotpath,cainfo=system.file("CurlSSL","cacert.pem",package="RCurl"))
                    tempGOdb<-strsplit(uniprothtml,"\n")[[1]]
                    for(EntryIDcount in 1:nrow(outputresult))
                       {
                          if(outputresult[EntryIDcount,"Uniprot Entry"]!="Not found" |outputresult[EntryIDcount,"Uniprot Entry"]!="More than 1")
                             {
                                 lookingpattern<-paste(outputresult[EntryIDcount,"Uniprot Entry"],"\t",sep="")
                                 eachtempGOdb<-grep(lookingpattern,tempGOdb,value=T)
                                 eachtempGOdb2<-gsub(lookingpattern,replacement="",eachtempGOdb)
                                 eachtempGOdb2<-union(eachtempGOdb2,eachtempGOdb2)
                                 eachtempGOdb3<-setdiff(eachtempGOdb2,c("GO:0003674","GO:0005575","GO:0008150"))
                                 if(length(eachtempGOdb3)!=0)
                                    {
                                       eachtotalGOID<-informationcombine(eachtempGOdb3,"; ")
                                       outputresult[EntryIDcount,"go-id"]<-eachtotalGOID
                                    }else
                                       {
                                           outputresult[EntryIDcount,"go-id"]<-""
                                       }
                             }
                       }
                }
             if(addintensityranking==TRUE)
                {
                   comparedata<-analyseIDresult[[filecount]][,c(1,which(colnames(analyseIDresult[[filecount]])=="Intensityranking"))]
                	 colnames(comparedata)<-c(searchcolumn2,"Intensityranking")
                   outputresult2<-merge(outputresult,comparedata,all.x=TRUE,by=searchcolumn2)
                   outputresult2<-outputresult2[,c(ncol(outputresult2),1:(ncol(outputresult2)-1))]
                   unionuniprotentry<-union(outputresult2[,"Uniprot Entry"],outputresult2[,"Uniprot Entry"])
                   for(x in 1:length(unionuniprotentry))
                      {
                          foundrearrangeposition<-which(unionuniprotentry[x]==outputresult2[,"Uniprot Entry"])
                          if(x==1)
                             {
                                 rearrangeposition<-foundrearrangeposition[1]
                             }else
                                {
                                    rearrangeposition<-c(rearrangeposition,foundrearrangeposition[1])
                                }
                      }
                   outputresult2<-outputresult2[rearrangeposition,]
                	 #UNKNOWN FUNCTION
                	 outputresult3<-outputresult2[order(as.numeric(outputresult2[,"Intensityranking"]),decreasing=FALSE),]
                	 analyseIDresult[[filecount]]<-outputresult3
                }else
                   {
                   	   analyseIDresult[[filecount]]<-outputresult
                   }

              writeDataTable(assayanalysisrsult, reportnames[filecount], analyseIDresult[[filecount]], startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = FALSE, tableStyle = "TableStyleLight9")

          }

       
       reportdirection<-gsub("-uniprot",replacement="",reportdirection)
       reportdirection<-gsub(".xlsx",replacement="",reportdirection)
       reportdirection<-paste(reportdirection,"-","uniprot",".xlsx",sep="")
       saveWorkbook(assayanalysisrsult, reportdirection, overwrite = TRUE)


   }


  # grep("YGL150C",organismidDBlist[[organismidDBcount]][1,"Gene.name"])
#============================================================================================================================================================================#
#collect protein seq(".fasta") from uniprot, this offer file to search common motif on MEME GLAM2

Uniprotseqsearch<-function(reportdirection)
   {
       report<-loadWorkbook(reportdirection)       
       assayanalysisrsult<-createWorkbook(creator =username())
       reportnames<-grep("hit",names(report),value=T)
       analyseIDresult<-vector("list",length(names(report)))
       for(filecount in 1:length(reportnames))
          {
             #讀取檔案，與製作新的表格
             addWorksheet(assayanalysisrsult, reportnames[filecount])
             analyseIDresult[[filecount]]<-readWorkbook(report, sheet = reportnames[filecount],startRow = 1,colNames = TRUE,rowNames = FALSE)
             for(getdatacounnt in 1:nrow(analyseIDresult[[filecount]]))
                {
                	uniprotentry<-as.character(analyseIDresult[[filecount]][getdatacounnt,"Uniprot.Entry"])
                	uniprotpath<-paste("http://www.uniprot.org/uniprot/",uniprotentry,".fasta",sep="")
                	uniprothtml <- getURL(uniprotpath,cainfo=system.file("CurlSSL","cacert.pem",package="RCurl"))
                	temp<-strsplit(uniprothtml,"\n")[[1]]
                	temp2<-c(temp,"")
                	result<-t(t(temp2))
                	if(getdatacounnt==1)
                	   {
                	   	output<-result
                	   }else
                	      {
                	      	 output<-rbind(output,result)
                	      }

                }
             output<-as.data.frame(output)
             writeData(assayanalysisrsult, reportnames[filecount], output, startCol = 1, startRow = 1, xy = NULL,colNames = FALSE,rowNames = FALSE)
          }
       reportdirection2<-gsub("-uniprot",replacement="",reportdirection)
       reportdirection2<-gsub(".xlsx",replacement="",reportdirection)
       reportdirection2<-paste(reportdirection,"-","uniprot sequence for GLAM2",".xlsx",sep="")
       saveWorkbook(assayanalysisrsult, reportdirection2, overwrite = TRUE)
   }


#============================================================================================================================================================================#

#three kind protein name can chose,"Systematic.name" ,"Gene.name", "Uniprot.Entry"

showupname<-"Gene.name"
GOsearch<-function(reportdirection,showupname)
   {

       report<-loadWorkbook(reportdirection)       
       assayanalysisrsult<-createWorkbook(creator =username())
       
       reportnames<-grep("hit",names(report),value=T)
       analyseIDresult<-vector("list",length(names(report)))
       for(filecount in 1:length(reportnames))
          {
             #讀取檔案，與製作新的表格
             addWorksheet(assayanalysisrsult, reportnames[filecount])
             analyseIDresult[[filecount]]<-readWorkbook(report, sheet = reportnames[filecount],startRow = 1,colNames = TRUE,rowNames = FALSE)
             goresultlist<-vector("list",length=nrow(analyseIDresult[[filecount]]))
             for(x in 1:nrow(analyseIDresult[[filecount]]))
                {
                	   goresultlist[[x]]<-strsplit(as.character(analyseIDresult[[filecount]][x,"go-id"]),"; ")[[1]]
                	   if(x==1)
                	      {
                	      	combine<-goresultlist[[x]]
                	      }else
                	         {
                	         	combine<-union(combine,goresultlist[[x]])

                	         }
                }
             combine<-combine[which(is.na(combine)!=TRUE)]


             combinecount<-combine
             outputresult<-select(GO.db, keys=combine, columns=c("TERM","ONTOLOGY","DEFINITION"),keytype="GOID")
             replicateGOIDproteinlist<-vector("list",length(combinecount))
             for(x in 1:length(combine))
                {
                	   combinecount[x]<-length(grep(combine[x],analyseIDresult[[filecount]][,"go-id"]))
                	   replicateGOIDproteinlist[[x]]<-analyseIDresult[[filecount]][grep(combine[x],analyseIDresult[[filecount]][,"go-id"]),showupname]
                }
             combinecount<-as.numeric(combinecount)
             replicateproteinmatrix<-matrix("",ncol=max(as.numeric(combinecount)),nrow=length(combine))

             for(x in 1:length(combine))
                {
                	   replicateproteinmatrix[x,1:length(replicateGOIDproteinlist[[x]])]<-replicateGOIDproteinlist[[x]]
                }
             searchlist2<-cbind(outputresult,replicateproteinmatrix)
             #summarycheck

             #aa<-read.csv("goresultwithinvolvedprotein.csv")
             #searchlist2<-aa[2:ncol(aa)]
             for(x in 1:ncol(searchlist2))
                {
                     searchlist2[,x]<-as.character(searchlist2[,x])
                }
             ontologykind<-c("BP","CC","MF")
             ontologyproteincount<-c(1:3)
             names(ontologyproteincount)<-ontologykind
             noinformationproteinlist<-vector("list",3)
             involvedproteinlist<-vector("list",3)
             summaryresult<-matrix("",ncol=6,nrow=3)
             colnames(summaryresult)<-c("Ontology","Ontology count","Involved protein count","Involved protein list","Not involved protein count","Not involved protein list")
             rownames(summaryresult)<-c(reportnames[filecount],"Totalhit",nrow(analyseIDresult[[filecount]]))

             for(ontologycount in 1:3)
                {
                  summaryresult[ontologycount,"Ontology"]<-ontologykind[ontologycount]
                  searchlist3<-searchlist2[which(searchlist2[,"ONTOLOGY"]==ontologykind[ontologycount]),5:14]
                  summaryresult[ontologycount,"Ontology count"]<-nrow(searchlist3)
                  for(x in 1:nrow(searchlist3))
                     {
                       #tempsave<-vector("numeric",length=length(which(searchlist3[x,]!="")))
                       #tempsave[1:length(which(searchlist3[x,]!=""))]<-searchlist3[x,which(searchlist3[x,]!="")]
                       if(x==1)
                          {
                            unionproteinlist<-as.character(searchlist3[x,])
                          }else
                             {
                                 unionproteinlist<-union(unionproteinlist,as.character(searchlist3[x,]))
                             }
                     }
                  unionproteinlist<-unionproteinlist[which(unionproteinlist!="")]
                  noinformationproteinlist[[ontologycount]]<-setdiff(as.character(analyseIDresult[[filecount]][,showupname]),unionproteinlist)
                  summaryresult[ontologycount,"Not involved protein list"]<-informationcombine(noinformationproteinlist[[ontologycount]],";")
                  summaryresult[ontologycount,"Not involved protein count"]<-length(noinformationproteinlist[[ontologycount]])
                  summaryresult[ontologycount,"Involved protein list"]<-informationcombine(unionproteinlist,";")
                  summaryresult[ontologycount,"Involved protein count"]<-length(unionproteinlist)
                }

             if(filecount==1)
                {
                  outputsummaryresult<-summaryresult
                }else
                   {
                       outputsummaryresult<-rbind(outputsummaryresult,summaryresult)
                   }
             #write.csv(searchlist2,paste(getwd(),"/","goresultwithinvolvedprotein.csv",sep=""))

             writeDataTable(assayanalysisrsult, reportnames[filecount], searchlist2, startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = FALSE, tableStyle = "TableStyleLight9")

          }
       outputsummaryresult<-as.data.frame(outputsummaryresult)
       addWorksheet(assayanalysisrsult, "Summary result")
       writeDataTable(assayanalysisrsult, "Summary result", outputsummaryresult, startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = TRUE, tableStyle = "TableStyleLight9")
       
       reportdirection2<-gsub("-uniprot",replacement="",reportdirection)
       reportdirection2<-gsub(".xlsx",replacement="",reportdirection2)
       reportdirection2<-paste(reportdirection2,"-","GO",".xlsx",sep="")
       saveWorkbook(assayanalysisrsult, reportdirection2, overwrite = TRUE)


   }



#============================================================================================================================================================================#
ancestorlayer<-3
#if searchedinputtarget is specific vector of GO ID it will search which children GOID(from protein) and protein is belong to it, and will close ancestorlayer function
#also,must give searchedinputtarget correct searchedinputtargetontologytype(ex:"BP","CC","MF").Remember, only can search 1 kind of ontology 1 time.
searchedinputtarget<-""
#searchedinputtargetontologytype is mean classify protein depend on previous GOID search result
searchedinputtargetontologytype<-"auto"

GOclassify<-function(reportdirection,uniprotreportdir,ancestorlayer,searchedinputtarget,searchedinputtargetontologytype)
   {

       GOreport<-loadWorkbook(GOreportdir)       
       #uniprotreport<-loadWorkbook(uniprotreportdir)

       assayanalysisrsult<-createWorkbook(creator =username())
       reportnames<-grep("hit",names(GOreport),value=T)
       reportnames<-grep("Summary result",reportnames,value=T,invert=TRUE)
       analyseIDresult<-vector("list",(length(names(GOreport))+1))
       #proteinlist<-analyseIDresult
       for(filecount in 1:length(reportnames))
          {
             #讀取檔案，與製作新的表格
               addWorksheet(assayanalysisrsult, reportnames[filecount])
               analyseIDresult[[filecount]]<-readWorkbook(GOreport, sheet = reportnames[filecount],startRow = 1,colNames = TRUE,rowNames = FALSE)

               #proteinlist[[filecount]]<-readWorkbook(uniprotreport, sheet = reportnames[filecount],startRow = 1,colNames = TRUE,rowNames = FALSE)[,c("Systematic.name" ,"Gene.name", "Uniprot.Entry")]
               analyseGOIDlist<-analyseIDresult[[filecount]][,"GOID"]
               inputGOproteinlist<-analyseIDresult[[filecount]][,5:ncol(analyseIDresult[[filecount]])]
               rownames(inputGOproteinlist)<-analyseGOIDlist
               #collect go data from GO.db
               # Convert the object to a list
               mainontologysearchvector<-c("BP","CC","MF")
             if(searchedinputtarget[1]=="")
                {
                	searchedinputtargetontologytype<-searchedinputtargetontologytype
                	searchedinputtargetontologytype<-"auto"
                	Start<-1
                	End<-3
                }else
                   {
                   	   Start<-which(mainontologysearchvector==searchedinputtargetontologytype)
                   	   End<-Start
                   }
             for(mainontologysearchcount in Start:End)
                {
                	if(mainontologysearchcount==1)
                	   {
                	   	   GOancestorlist <- as.list(GOBPANCESTOR)
                	   }
                	if(mainontologysearchcount==2)
                	   {
                	   	   GOancestorlist <- as.list(GOCCANCESTOR)
                	   }
                	if(mainontologysearchcount==3)
                	   {
                	   	   GOancestorlist <- as.list(GOMFANCESTOR)
                	   }

                	# Remove GO IDs that do not have any ancestor
                	GOancestorlist <- GOancestorlist[!is.na(GOancestorlist)]
                	GOancestornamelist<-names(GOancestorlist)
                	GOancestorlengthcount<-unlist(lapply(GOancestorlist,length))-1
                	#資料抓取時，取母層只有一層或兩層的做分類即可，先選用一層地做分類
                	if(searchedinputtarget[1]=="")
                	   {
                	   	   if(ancestorlayer!=1)
                            {
                                
                                for(ancestorlayer2 in 1:ancestorlayer)
                                   {
                                       totalclassposition<-which(GOancestorlengthcount==ancestorlayer2)
                                       if(ancestorlayer2==1)
                                          {
                                              totalclassposition2<-totalclassposition
                                          }else
                                             {
                                                 totalclasspositionacstor<-vector("list",length(totalclassposition))
                                                 for(totalclasspositionacstorcount in 1:length(totalclassposition))
                                                    {
                                                        totalclasspositionacstor[[totalclasspositionacstorcount]]<-GOancestorlist[[totalclassposition[totalclasspositionacstorcount]]]
                                                        replicateresult<-intersect(totalclasspositionacstor[[totalclasspositionacstorcount]],names(totalclassposition2))
                                                        if(totalclasspositionacstorcount==1)
                                                           {
                                                              replicateresult2<-replicateresult
                                                           }else
                                                              {
                                                                  replicateresult2<-union(replicateresult2,replicateresult)
                                                              }
                                                    }
                                                totalclassposition2<-c(totalclassposition2[setdiff(names(totalclassposition2),replicateresult2)],totalclassposition)
                                                totalclassposition2<-totalclassposition2[order(as.numeric(totalclassposition2),decreasing=FALSE)]


                                             }
                                   }
                                totalclassposition<-totalclassposition2
                            }else
                               {
                                   totalclassposition<-which(GOancestorlengthcount==ancestorlayer)
                               }
                         
                	   }else
                	      {
                	          totalclassposition<-vector("numeric",length(searchedinputtarget))
                	          names(totalclassposition)<-searchedinputtarget
                	          for(searchedinputtargetcount in 1:length(searchedinputtarget))
                	             {
                	             	 totalclassposition[searchedinputtargetcount]<-which(GOancestornamelist==searchedinputtarget[searchedinputtargetcount])
                	             }
                	          totalclassposition<-totalclassposition[order(as.numeric(totalclassposition),decreasing=FALSE)]
                	          
                	      }
                	
                	analyseseprateontologytarget<-intersect(analyseGOIDlist,names(GOancestorlist))
                	analyseseprateontologytargetlist<-vector("list",length=length(analyseseprateontologytarget))
                	totalclassproteinlist<-vector("list",length=length(totalclassposition))
                	names(totalclassproteinlist)<-names(totalclassposition)
                	totalclassGOIDlist<-totalclassproteinlist
                	#比對ancestor，將有分類目標的部份擷取出來
                	for(searchcount in 1:length(analyseseprateontologytargetlist))
                	   {
                	   	searchingGOID<-analyseseprateontologytarget[searchcount]
                	   	foundresult<-intersect(names(totalclassposition),GOancestorlist[[searchingGOID]])

                	   	#foundresult 是 最後分類的目標target，不是所搜尋的GOID
                	   	if(length(foundresult)!=0)
                	   	   {
                	   	   	   analyseseprateontologytargetlist[[searchcount]]<-as.character(foundresult)
                	   	   	   for(foundresultcount in 1:length(foundresult))
                	   	   	      {
                	   	   	      	  searchedclassproteinlist<-as.character(inputGOproteinlist[searchingGOID,])
                	   	   	      	  searchedclassproteinlist<-searchedclassproteinlist[which(searchedclassproteinlist!="")]
                	   	   	      	  outputtargetclass<-foundresult[foundresultcount]
                	   	   	      	  if(length(totalclassproteinlist[[outputtargetclass]]!=0))
                	   	   	      	     {
                	   	   	      	         previousclassplrteinlistresult<-totalclassproteinlist[[outputtargetclass]]
                	   	   	      	         totalclassproteinlist[[outputtargetclass]]<-union(previousclassplrteinlistresult,searchedclassproteinlist)

                	   	   	      	         previousclassGOIDlistresult<-totalclassGOIDlist[[outputtargetclass]]
                	   	   	      	         totalclassGOIDlist[[outputtargetclass]]<-union(previousclassGOIDlistresult,searchingGOID)
                	   	   	      	     }else
                	   	   	      	        {
                	   	   	      	        	totalclassproteinlist[[outputtargetclass]]<-searchedclassproteinlist

                	   	   	      	        	totalclassGOIDlist[[outputtargetclass]]<-searchingGOID
                	   	   	      	        }
                	   	   	      }
                	   	   }
                	   }
                	totalclassinformation<-select(GO.db, keys=names(totalclassposition), columns=c("TERM","ONTOLOGY","DEFINITION"),keytype="GOID")
                	totalclassproteinlistcount<-unlist(lapply(totalclassproteinlist,length))
                	splitunit<-";"
                	combineclassproteinlist<-informationlistcombine(totalclassproteinlist,splitunit)

                	totalclassGOIDlistcount<-unlist(lapply(totalclassGOIDlist,length))
                	combineclassGOIDlist<-informationlistcombine(totalclassGOIDlist,splitunit)
                	summarizedresult<-cbind(totalclassproteinlistcount,combineclassproteinlist,totalclassGOIDlistcount,combineclassGOIDlist)
                	colnames(summarizedresult)<-c("Protein count","Protein list","GOID count","GOID list")
                	outputresult<-cbind(totalclassinformation,summarizedresult)
                	outputresult2<-outputresult[which(outputresult[,"Protein count"]!=0),]
                	if(mainontologysearchcount==1)
                	   {
                	   	   finaloutput<-outputresult2
                	   }else
                	      {
                	      	  finaloutput<-rbind(finaloutput,outputresult2)
                	      }


                }
             
             writeDataTable(assayanalysisrsult, reportnames[filecount], finaloutput, startCol = 1, startRow = 1, xy = NULL,colNames = TRUE,rowNames = FALSE, tableStyle = "TableStyleLight9")

          }

       
       reportdirection2<-gsub("-GO",replacement="",GOreportdir)
       reportdirection2<-gsub(".xlsx",replacement="",reportdirection2)
       if(searchedinputtarget[1]=="")
          {
              reportdirection2<-paste(GOreportdir,"-","GO classification"," with ",ancestorlayer," ancestor layer",".xlsx",sep="")
          }else
             {
                 reportdirection2<-paste(GOreportdir,"-","GO classification"," with ","searchedinputtarget",".xlsx",sep="")
             }
       
       saveWorkbook(assayanalysisrsult, reportdirection2, overwrite = TRUE)


   }
#============================================================================================================================================================================#
