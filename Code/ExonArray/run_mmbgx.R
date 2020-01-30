rm(list=ls());
setwd("/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/ExonArray/Data3");
library(oligo);
library(mmbgx); # load MMGBX library
library(affyio);

d <- read.celfiles(dir(),pkgname="pd.huex.1.0.st.v2"); # run mmbgx to generate isoform expression
mmbgx(d, arrayType="HuEx-1_0-st-v2",geneLevel=FALSE, inputDirs=paste("transcriptInput", c(1:2), sep="."), standalone=TRUE);
# run mmbgx to generate gene expression
mmbgx(d, arrayType="HuEx-1_0-st-v2", geneLevel=TRUE, inputDirs=paste("geneInput", c(1:2), sep="."), standalone=TRUE);

#combineRuns.mmbgx(dir("geneRuns1", full.names=TRUE), sampleSets=1, outputDir="gene1");
resG1 <- readSingle.mmbgx("geneRuns1/run.1");
X1 = cbind(resG1$muave,resG1$PSids);
rm(resG1);

resG2 <- readSingle.mmbgx("geneRuns2/run.2");
X2 = cbind(resG2$muave,resG2$PSids);
rm(resG2);

resG3 <- readSingle.mmbgx("geneRuns3/run.3");
X3 = cbind(resG3$muave,resG3$PSids);
rm(resG3);

resG4 <- readSingle.mmbgx("geneRuns4/run.4");
X4 = cbind(resG4$muave,resG4$PSids);
rm(resG4);

resG5 <- readSingle.mmbgx("geneRuns5/run.5");
X5 = cbind(resG5$muave,resG5$PSids);
rm(resG5);

resG6 <- readSingle.mmbgx("geneRuns6/run.6");
X6 = cbind(resG6$muave,resG6$PSids);
rm(resG6);

resG7 <- readSingle.mmbgx("geneRuns7/run.7");
X7 = cbind(resG7$muave,resG7$PSids);
rm(resG7);

resG8 <- readSingle.mmbgx("geneRuns8/run.8");
X8 = cbind(resG8$muave,resG8$PSids);
rm(resG8);

Data = cbind(X1[,1],X2[,1],X3[,1],X4[,1],X5[,1],X6[,1],X7[,1],X8[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('BT549','HCC1937','Hs578T','MCF7','MDA-MB-231','MDA-MB-436','SK-BR-3','T47D');
write.table(Data, file = "Gene1.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");

resT1 <- readSingle.mmbgx("transcriptRuns1/run.1");
X1 = cbind(resT1$muave,resT1$PSids);
G = resT1$geneIDs;
G_T = resT1$geneTranscripts;
rm(resT1);

resT2 <- readSingle.mmbgx("transcriptRuns2/run.2");
X2 = cbind(resT2$muave,resT2$PSids);
rm(resT2);

resT3 <- readSingle.mmbgx("transcriptRuns3/run.3");
X3 = cbind(resT3$muave,resT3$PSids);
rm(resT3);

resT4 <- readSingle.mmbgx("transcriptRuns4/run.4");
X4 = cbind(resT4$muave,resT4$PSids);
rm(resT4);

resT5 <- readSingle.mmbgx("transcriptRuns5/run.5");
X5 = cbind(resT5$muave,resT5$PSids);
rm(resT5);

resT6 <- readSingle.mmbgx("transcriptRuns6/run.6");
X6 = cbind(resT6$muave,resT6$PSids);
rm(resT6);

resT7 <- readSingle.mmbgx("transcriptRuns7/run.7");
X7 = cbind(resT7$muave,resT7$PSids);
rm(resT7);

resT8 <- readSingle.mmbgx("transcriptRuns8/run.8");
X8 = cbind(resT8$muave,resT8$PSids);
rm(resT8);


Data = cbind(X1[,1],X2[,1],X3[,1],X4[,1],X5[,1],X6[,1],X7[,1],X8[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('BT549','HCC1937','Hs578T','MCF7','MDA-MB-231','MDA-MB-436','SK-BR-3','T47D');
write.table(Data, file = "Transcript1.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");


setwd("/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/ExonArray/Data2");

resT1 <- readSingle.mmbgx("TranscriptRuns1/run.1");
X1 = cbind(resT1$muave,resT1$PSids);
G = resT1$geneIDs;
G_T = resT1$geneTranscripts;
rm(resT1);

resT2 <- readSingle.mmbgx("TranscriptRuns2/run.2");
X2 = cbind(resT2$muave,resT2$PSids);
rm(resT2);

resT3 <- readSingle.mmbgx("TranscriptRuns3/run.3");
X3 = cbind(resT3$muave,resT3$PSids);
rm(resT3);

resT4 <- readSingle.mmbgx("TranscriptRuns4/run.4");
X4 = cbind(resT4$muave,resT4$PSids);
rm(resT4);

resT5 <- readSingle.mmbgx("TranscriptRuns5/run.5");
X5 = cbind(resT5$muave,resT5$PSids);
rm(resT5);

resT6 <- readSingle.mmbgx("TranscriptRuns6/run.6");
X6 = cbind(resT6$muave,resT6$PSids);
rm(resT6);

resT7 <- readSingle.mmbgx("TranscriptRuns7/run.7");
X7 = cbind(resT7$muave,resT7$PSids);
rm(resT7);

resT8 <- readSingle.mmbgx("TranscriptRuns8/run.8");
X8 = cbind(resT8$muave,resT8$PSids);
rm(resT8);

resT9 <- readSingle.mmbgx("TranscriptRuns9/run.9");
X9 = cbind(resT9$muave,resT9$PSids);
rm(resT9);

resT10 <- readSingle.mmbgx("TranscriptRuns10/run.10");
X10 = cbind(resT10$muave,resT10$PSids);
rm(resT10);

resT11 <- readSingle.mmbgx("TranscriptRuns11/run.11");
X11 = cbind(resT11$muave,resT11$PSids);
rm(resT11);

resT12 <- readSingle.mmbgx("TranscriptRuns12/run.12");
X12 = cbind(resT12$muave,resT12$PSids);
rm(resT12);

resT13 <- readSingle.mmbgx("TranscriptRuns13/run.13");
X13 = cbind(resT13$muave,resT13$PSids);
rm(resT13);

resT14 <- readSingle.mmbgx("TranscriptRuns14/run.14");
X14 = cbind(resT14$muave,resT14$PSids);
rm(resT14);

resT15 <- readSingle.mmbgx("TranscriptRuns15/run.15");
X15 = cbind(resT15$muave,resT15$PSids);
rm(resT15);

resT16 <- readSingle.mmbgx("TranscriptRuns16/run.16");
X16 = cbind(resT16$muave,resT16$PSids);
rm(resT16);

resT17 <- readSingle.mmbgx("TranscriptRuns17/run.17");
X17 = cbind(resT17$muave,resT17$PSids);
rm(resT17);

resT18 <- readSingle.mmbgx("TranscriptRuns18/run.18");
X18 = cbind(resT18$muave,resT18$PSids);
rm(resT18);

resT19 <- readSingle.mmbgx("TranscriptRuns19/run.19");
X19 = cbind(resT19$muave,resT19$PSids);
rm(resT19);

Data = cbind(X1[,1],X2[,1],X3[,1],X4[,1],X5[,1],X6[,1],X7[,1],X8[,1],X9[,1],X10[,1],X11[,1],X12[,1],X13[,1],X14[,1],X15[,1],X16[,1],X17[,1],X18[,1],X19[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('PANC-1','H460','OVCAR-3','OVCAR-4','OVCAR-8','SK-OV-3','DU-145','PC-3','Hela','AGS','HT1080','IMR90','WI38','HCT-116','HCT-15','HT29','KM12','SW-620','A549');
write.table(Data, file = "Transcript2.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");


setwd("/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/ExonArray/Data2");

resT1 <- readSingle.mmbgx("geneRuns1/run.1");
X1 = cbind(resT1$muave,resT1$PSids);
rm(resT1);

resT2 <- readSingle.mmbgx("geneRuns2/run.2");
X2 = cbind(resT2$muave,resT2$PSids);
rm(resT2);

resT3 <- readSingle.mmbgx("geneRuns3/run.3");
X3 = cbind(resT3$muave,resT3$PSids);
rm(resT3);

resT4 <- readSingle.mmbgx("geneRuns4/run.4");
X4 = cbind(resT4$muave,resT4$PSids);
rm(resT4);

resT5 <- readSingle.mmbgx("geneRuns5/run.5");
X5 = cbind(resT5$muave,resT5$PSids);
rm(resT5);

resT6 <- readSingle.mmbgx("geneRuns6/run.6");
X6 = cbind(resT6$muave,resT6$PSids);
rm(resT6);

resT7 <- readSingle.mmbgx("geneRuns7/run.7");
X7 = cbind(resT7$muave,resT7$PSids);
rm(resT7);

resT8 <- readSingle.mmbgx("geneRuns8/run.8");
X8 = cbind(resT8$muave,resT8$PSids);
rm(resT8);

resT9 <- readSingle.mmbgx("geneRuns9/run.9");
X9 = cbind(resT9$muave,resT9$PSids);
rm(resT9);

resT10 <- readSingle.mmbgx("geneRuns10/run.10");
X10 = cbind(resT10$muave,resT10$PSids);
rm(resT10);

resT11 <- readSingle.mmbgx("geneRuns11/run.11");
X11 = cbind(resT11$muave,resT11$PSids);
rm(resT11);

resT12 <- readSingle.mmbgx("geneRuns12/run.12");
X12 = cbind(resT12$muave,resT12$PSids);
rm(resT12);

resT13 <- readSingle.mmbgx("geneRuns13/run.13");
X13 = cbind(resT13$muave,resT13$PSids);
rm(resT13);

resT14 <- readSingle.mmbgx("geneRuns14/run.14");
X14 = cbind(resT14$muave,resT14$PSids);
rm(resT14);

resT15 <- readSingle.mmbgx("geneRuns15/run.15");
X15 = cbind(resT15$muave,resT15$PSids);
rm(resT15);

resT16 <- readSingle.mmbgx("geneRuns16/run.16");
X16 = cbind(resT16$muave,resT16$PSids);
rm(resT16);

resT17 <- readSingle.mmbgx("geneRuns17/run.17");
X17 = cbind(resT17$muave,resT17$PSids);
rm(resT17);

resT18 <- readSingle.mmbgx("geneRuns18/run.18");
X18 = cbind(resT18$muave,resT18$PSids);
rm(resT18);

resT19 <- readSingle.mmbgx("geneRuns19/run.19");
X19 = cbind(resT19$muave,resT19$PSids);
rm(resT19);

Data = cbind(X1[,1],X2[,1],X3[,1],X4[,1],X5[,1],X6[,1],X7[,1],X8[,1],X9[,1],X10[,1],X11[,1],X12[,1],X13[,1],X14[,1],X15[,1],X16[,1],X17[,1],X18[,1],X19[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('PANC-1','H460','OVCAR-3','OVCAR-4','OVCAR-8','SK-OV-3','DU-145','PC-3','Hela','AGS','HT1080','IMR90','WI38','HCT-116','HCT-15','HT29','KM12','SW-620','A549');
write.table(Data, file = "Gene2.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");

setwd("/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/ExonArray/Data3");

resG1 <- readSingle.mmbgx("geneRuns1/run.1");
X1 = cbind(resG1$muave,resG1$PSids);
rm(resG1);

resG2 <- readSingle.mmbgx("geneRuns2/run.2");
X2 = cbind(resG2$muave,resG2$PSids);
rm(resG2);

resG3 <- readSingle.mmbgx("geneRuns3/run.3");
X3 = cbind(resG3$muave,resG3$PSids);
rm(resG3);

resG4 <- readSingle.mmbgx("geneRuns4/run.4");
X4 = cbind(resG4$muave,resG4$PSids);
rm(resG4);

resG5 <- readSingle.mmbgx("geneRuns5/run.5");
X5 = cbind(resG5$muave,resG5$PSids);
rm(resG5);

resG6 <- readSingle.mmbgx("geneRuns6/run.6");
X6 = cbind(resG6$muave,resG6$PSids);
rm(resG6);

Data = cbind(X1[,1],X2[,1],X3[,1],X4[,1],X5[,1],X6[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('MCF10A','DLD1','A2780','Caov-3','ES-2','TOV-21G');
write.table(Data, file = "Gene3.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");



setwd("/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/ExonArray/Data3");

resT1 <- readSingle.mmbgx("TranscriptRuns1/run.1");
X1 = cbind(resT1$muave,resT1$PSids);
G = resT1$geneIDs;
G_T = resT1$geneTranscripts;
rm(resT1);

resT2 <- readSingle.mmbgx("TranscriptRuns2/run.2");
X2 = cbind(resT2$muave,resT2$PSids);
rm(resT2);

resT3 <- readSingle.mmbgx("TranscriptRuns3/run.3");
X3 = cbind(resT3$muave,resT3$PSids);
rm(resT3);

resT4 <- readSingle.mmbgx("TranscriptRuns4/run.4");
X4 = cbind(resT4$muave,resT4$PSids);
rm(resT4);

resT5 <- readSingle.mmbgx("TranscriptRuns5/run.5");
X5 = cbind(resT5$muave,resT5$PSids);
rm(resT5);

resT6 <- readSingle.mmbgx("TranscriptRuns6/run.6");
X6 = cbind(resT6$muave,resT6$PSids);
rm(resT6);

Data = cbind(X1[,1],X2[,1],X3[,1],X4[,1],X5[,1],X6[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('MCF10A','DLD1','A2780','Caov-3','ES-2','TOV-21G');
write.table(Data, file = "Transcript3.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");




setwd("/project/compbioRAID1/WeiZhang/ProbeDesign/NanostringData/ExonArray/Data4");

resG1 <- readSingle.mmbgx("geneRuns1/run.1");
X1 = cbind(resG1$muave,resG1$PSids);
rm(resG1);

resG2 <- readSingle.mmbgx("geneRuns2/run.2");
X2 = cbind(resG2$muave,resG2$PSids);
rm(resG2);

Data = cbind(X1[,1],X2[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('Hep_G2','Caco2');
write.table(Data, file = "Gene4.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");



resT1 <- readSingle.mmbgx("transcriptRuns1/run.1");
X1 = cbind(resT1$muave,resT1$PSids);
G = resT1$geneIDs;
G_T = resT1$geneTranscripts;
rm(resT1);

resT2 <- readSingle.mmbgx("transcriptRuns2/run.2");
X2 = cbind(resT2$muave,resT2$PSids);
rm(resT2);

Data = cbind(X1[,1],X2[,1]);
rownames(Data) = X1[,2];
colnames(Data) = c('Hep_G2','Caco2');
write.table(Data, file = "Transcript4.txt", append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "");
