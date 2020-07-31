setwd("/home/aurore/Documents/Rainette_FGB/VORTEX10/VOutput")
setwd("~/Aurore/VORTEX10")
#### Exemple of Campbell and Pacioni data ####
require(vortexRdata)
data("pac.clas")
head(pac.clas)
#raw data guaranty to work
pac.dir <- system.file("extdata", "pacioni", package="vortexRdata")
cam.dir <- system.file("extdata", "campbell", package="vortexRdata")
f <- file.path(pac.dir, "Pacioni_et_al_ST_Classic(Base).stdat")
require(vortexR)
one.st.classic <- collate_one_dat(f, 3)
head(one.st.classic)
data("pac.clas.Nadults")
head(pac.clas.Nadults)
# Function Dot plots of mean Vortex parameters
data(pac.clas)
dot <- dot_plot(data=pac.clas, project='Pacioni_et_al', scenario='ST_Classic',
                yrs=c(80, 120),
                params=c('PExtinct', 'Nextant', 'Het', 'Nalleles'),
                save2disk=FALSE)
# Collate Vortex .run output files
run <- collate_run('Pacioni_et_al', 'ST_LHS', 1, dir_in=pac.dir,
                   save2disk=FALSE)
# lookup_table Summary table of simulation parameters (pour les analyses de sensibilite seulement)
data(pac.clas)
lkup.st.classic <- lookup_table(data=pac.clas, project='Pacioni_et_al',
scenario='ST_Classic', pop='Population 1',SVs=c('SV1', 'SV2', 'SV3', 'SV4', 'SV5', 'SV6', 'SV7'),save2disk=FALSE)
# Funtion collate_yr before Nadults

# Nadults.
data(pac.yr)
NadultAll <- Nadults(data=pac.yr[[2]], scenarios='all', gen=2.54, yr0=50,
                     yrt=120, save2disk=FALSE)
------------------------------
#### Packages ####
install.packages("vortexR", dependencies = TRUE)
install.packages("devtools")
devtools::install_github("carlopacioni/vortexRdata")
devtools::install_github("carlopacioni/vortexR")
library(vortexR) #https://rdrr.io/cran/vortexR/
require(vortexRdata)
require(dplyr) #gestion data.frame : https://statistique-et-logiciel-r.com/initiation-a-la-manipulation-de-donnees-avec-le-package-dplyr/
require(ggplot2)
require(glmulti)
require(MASS) # when overdispersion ;negative binomial error distribution
require(pscl) # zero inflated model to deal with the excess of zeros caused by several parameter combinations
require(forcats) # modifier l ordre d un plot
#### Astuces eviter conflict packages ####
library(conflicted)
conflict_prefer("filter", "dplyr")
## [conflicted] Will prefer dplyr::filter over any other package

#### DATA extractions ####
# 1. Collate_DAT 
#Ajout des Basic scenario avec deux populations (+reintroduite)
basic_nondd <- collate_one_dat('Parameters-tests-nondd-RFGB_Basic nondd.dat', 1000, verbose = FALSE)
basic_dd <- collate_one_dat('Parameters-tests-RFGB_Basic dd.dat', 1000, verbose = FALSE)
allscenario_dd <- collate_dat('Parameters-tests-RFGB',1000, dir_in = "VOutput")
scenarioS89 <- collate_dat('Parameters-tests-RFGB_S89',1000, dec_sep=",", dir_in = "VOutput")
names(scenarioS89)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")


allscenario_nondd <- collate_dat('Parameters-tests-nondd-RFGB',1000, dir_in = "VOutput")

#Prévenir des bugs, renommer colomnes
names(basic_nondd)[9] <- "r.stoch"
names(basic_dd)[9] <- "r.stoch"
save(basic_dd,file="Basic-dd-2pops.rda")
save(basic_nondd,file="Basic-nondd-2pops.rda")

allscenario_dd <- read.table("VortexR/allscenario_dd - Copie.txt",sep=';') #sans Basicdd ajout après
merge-dd<-union(allscenario_dd,basic_dd,by.y=38)#trouver solution pour que les 2data.frame soit assemblée avant rRec pour reintroduite
allscenario_nondd <- read.table("VortexR/allscenario_nondd.txt",sep=';')
merge-nondd<-union(allscenario_nondd,basic_nondd)

# Commparaison scenarios reintroduction DD
reintro.3ans <- collate_dat('catastrophy-dd_3',1000, dec_sep=",", dir_in = "VOutput")
names(reintro.3ans)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
names(reintro.3ans)[9] <- "r.stoch"
reintro.3ans <- reintro.3ans %>% filter(pop.name %in% "Reintroduite") 

reintro.5ans <- collate_dat('catastrophy-dd_5',1000, dec_sep=",", dir_in = "VOutput")
names(reintro.5ans)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
names(reintro.5ans)[9] <- "r.stoch"
reintro.5ans <- reintro.5ans %>% filter(pop.name %in% "Reintroduite") 

# Catastrpohy
#avec meilleurs combinaisons
setwd("~/Aurore/VORTEX10/VOutput")
basic.cata25 <- collate_dat('catastrophy-dd_Basic-25ans', 1000, verbose = FALSE, dec_sep=",")
names(basic.cata25)[9] <- "r.stoch"
names(basic.cata25)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
cata25 <- collate_dat('catastrophy-dd_25ans',1000, verbose = FALSE, dec_sep=",")
names(cata25)[9] <- "r.stoch"
names(cata25)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
basic.cata50 <- collate_dat('catastrophy-dd_Basic-50ans', 1000, verbose = FALSE, dec_sep=",")
names(basic.cata50)[9] <- "r.stoch"
names(basic.cata50)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
cata50 <- collate_dat('catastrophy-dd_50ans',1000, verbose = FALSE, dec_sep=",")
names(cata50)[9] <- "r.stoch"
names(cata50)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")

setwd("~/Aurore/VORTEX10/VortexR")
save(allcata,file="allcata.rda")

#avec pires combinaisons
basic.worst25 <- collate_dat('catastrophy-worst_Basic-25ans', 1000, verbose = FALSE, dec_sep=",")
names(basic.worst25)[9] <- "r.stoch"
names(basic.worst25)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
cata.worst25 <- collate_dat('catastrophy-worst_25ans',1000, verbose = FALSE, dec_sep=",")
names(cata.worst25)[9] <- "r.stoch"
names(cata.worst25)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
basic.worst50 <- collate_dat('catastrophy-worst_Basic-50ans', 1000, verbose = FALSE, dec_sep=",")
names(basic.worst50)[9] <- "r.stoch"
names(basic.worst50)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")
cata.worst50 <- collate_dat('catastrophy-worst_50ans',1000, verbose = FALSE, dec_sep=",")
names(cata.worst50)[9] <- "r.stoch"
names(cata.worst50)[24:36]<-c("Nalleles","SE.Nalleles.","SD.Nalleles.","Nlethals","SE.Nlethals.","SD.Nlethals.","TE","SE.TE.","SD.TE.","medTE","MortalityT","SE.MortalityT.","SD.MortalityT.")

# Avec basic scenario a une seule population = naturelle (original version)
allscenario_dd <- read.table("allscenario_dd.txt",sep=';') # avec Basic 
allscenario_nondd <- read.table("allscenario_nondd.txt",sep=';') # avec Basic 

#Prévenir des bugs, renommer colomnes
names(allscenario_dd)[1] <- "scen.name"
names(allscenario_nondd)[1] <- "scen.name"
names(allscenario_dd)[9] <- "r.stoch"
names(allscenario_nondd)[9] <- "r.stoch"

# 2. Collate_run
run.dd <- collate_run('Parameters-tests-RFGB',scenario='', npops=3, dir_in = "VOutput", dec_sep=',')
basicrun.dd <- collate_run('Parameters-tests-RFGB',scenario='Basic', npops=3, dir_in = "VOutput", dec_sep=',')
run.nondd <- collate_run('Parameters-tests-nondd-RFGB',scenario='', npops=3, dir_in = "VOutput", dec_sep=',')
basicrun.nondd <- collate_run('Parameters-tests-RFGB',scenario='Basic', npops=3, dir_in = "VOutput", dec_sep=',')

setwd("~/Aurore/VORTEX10/catastrophy-dd")
run.cata <- collate_run('catastrophy-dd',scenario='', npops=3, dir_in = "VOutput", dec_sep=',')

# 3. Tables des valeurs moyennes des scenarios / a changer pour as.numeric
sum.dd <- read.csv("sorties-tests-parametres.csv")
sum.nondd <- read.csv("sorties-tests-parametres-nondd.csv")
# Split des populations
split.sum.dd <- split(sum.dd, sum.dd$Population)
split.sum.nondd <- split(sum.nondd, sum.nondd$Population)


#### Plot DATA ####
# ATTENTION ce ne sont pas des moyennes cumulées mais bien les valeurs de l annee plottee simplement
rec.dd <- dot_plot(data=best.dd.nat, project="Parameters-tests-RFGB.xlm", scenario='Basic dd',yrs=c(7,20),
               params=c("r.stoch","Nall","Het"), setcolour="Nall",plotpops='Naturelle',save2disk=FALSE)
rec.dd2 <- dot_plot(data=best.dd.reintro, project="Parameters-tests-RFGB.xlm", scenario='S',yrs=c(7,20),
                   params=c("r.stoch","Nall","Het"), setcolour="Nall",plotpops='Reintroduite',save2disk=FALSE)

rec.nondd1 <- dot_plot(data=allscenario_nondd, project="/home/aurore/Documents/Rainette_FGB/VORTEX10/Parameters-tests-nondd-RFGB.xlm", scenario='Basic dd',
                   yrs=c(10, 25), params="Nall", plotpops='Naturelle',save2disk=FALSE)
# Plot Pextant
names(allscenario_dd)[37] <- "SD.Pextant."
names(allscenario_dd)[38] <- "SD.Pextinct."

# Toutes les valeurs avec SD --» remarque des clusters et separations entre scenarios
dot_plot(data=allscenario_dd, project="Parameters-tests-RFGB.xlm", scenario='Basic dd',yrs=c(10,20),
         params="Pextant", setcolour="Nall",plotpops='Naturelle',save2disk=FALSE)
dot_plot(data=allscenario_dd, project="Parameters-tests-RFGB.xlm", scenario='Basic dd',yrs=c(10,20),
         params="Pextant", setcolour="Nall",plotpops='Reintroduite',save2disk=FALSE)
line_plot_year(data=Het.reintro.nondd, project="Parameters-tests-RFGB.xlm", scenario="S", params="GeneDiv", plotpops = "Reintroduite",save2disk = F)

#### Filtrer meilleurs en fonction de Pextant et Het les scenarios #####
Pextant.dd <- allscenario_dd %>% filter(Pextant>0.8)
Pextant.S89 <- scenarioS89 %>% filter(PExtant>0.8)

Het.natdd <- Pextant.dd %>% filter(pop.name %in% "Naturelle") %>%filter(Het>0.85)
Het.S89 <- scenarioS89 %>% filter(pop.name %in% "Naturelle") %>%filter(GeneDiv>0.85)
Het.reintrodd <- Pextant.dd %>% filter(pop.name %in% "Reintroduite") %>%filter(Het>0.90)

Pextant.nondd <- allscenario_nondd %>% filter(PExtant>0.8)
Het.nat.nondd <- Pextant.nondd %>% filter(pop.name %in% "Naturelle") %>%filter(GeneDiv>0.90)
Het.reintro.nondd <- Pextant.nondd %>% filter(pop.name %in% "Reintroduite") %>%filter(GeneDiv>0.90)
names(Het.nat.nondd)[20] <- "SD.GeneDiv."
names(Het.reintro.nondd)[20] <- "SD.GeneDiv."

dot_plot(data=Het.nat.nondd, project="Parameters-tests-nondd-RFGB.xlm", scenario='Basic nondd',yrs=c(10,25),
  params="GeneDiv", setcolour="Nextant", plotpops='Naturelle',save2disk=FALSE)
dot_plot(data=Pextant, project="Parameters-tests-RFGB.xlm", scenario='Basic dd',yrs=c(10,20),
         params="Pextant", plotpops='Reintroduite',save2disk=FALSE)

#### Filtrer pires en fonction de Pextant et Het les scenarios #####
Pextant.dd.worst <- allscenario_dd %>% filter(Pextant<0.3)
Het.natdd.worst <- Pextant.dd.worst %>% filter(pop.name %in% "Naturelle") %>%filter(Het<0.6)
tapply(Het.natdd.worst$Year, Het.natdd.worst$scen.name, length)
Het.reintrodd.worst <- Pextant.dd.worst %>% filter(pop.name %in% "Reintroduite") %>%filter(Het<0.60)
tapply(Het.reintrodd.worst$Year, Het.reintrodd.worst$scen.name, length)

dot_plot(data=Het.nat.nondd, project="Parameters-tests-nondd-RFGB.xlm", scenario='Basic nondd',yrs=c(10,25),
         params="GeneDiv", setcolour="Nextant", plotpops='Naturelle',save2disk=FALSE)
dot_plot(data=Pextant, project="Parameters-tests-RFGB.xlm", scenario='Basic dd',yrs=c(10,20),
         params="Pextant", plotpops='Reintroduite',save2disk=FALSE)

#### Comparer combinaisons DD reintro et naturelle par distri Het ####
best.Het.5ans <- Het.natdd %>% filter(scen.name %in% c("S1","S17","S9","S25","Basic dd"))
best.Het.3ans <- Het.natdd %>% filter(scen.name %in% c("S18","S2","S10","S26","Basic dd"))
best.Het.813 <- Het.natdd %>% filter(scen.name %in% c("S9","S1","S17","S25","Basic dd"))
best.Het.2040 <- Het.natdd %>% filter(scen.name %in% c("S89","S83","S81","S75","Basic dd"))

mean.het3ans <- na.omit(as.data.frame(tapply(best.Het.3ans$Nextant, best.Het.3ans$scen.name,mean)))
names(mean.het3ans)[1] <- "mean"
mean.het5ans <- na.omit(as.data.frame(tapply(best.Het.5ans$Nextant, best.Het.5ans$scen.name,mean)))
names(mean.het5ans)[1] <- "mean"
ks.test(mean.het3ans$mean,mean.het5ans$mean, alternative="greater") #p-value = 0.01832 donc 5 ans est bien supérieur à 3 ans

# Probleme avec cette methode car on calcul les moyennes mais pas les ecarts types et si ils se chevauchent finalement les tests 
# ne comptent plus... a verifier

mean.het813 <- na.omit(as.data.frame(tapply(best.Het.813$Pextant, best.Het.813$scen.name,mean)))
names(mean.het813)[1] <- "mean"
mean.het2040 <- na.omit(as.data.frame(tapply(Het.S89$PExtant, Het.S89$scen.name,mean)))
names(mean.het2040)[1] <- "mean"
ks.test(mean.het813$mean,mean.het2040$mean, alternative="two.sided") #p-value = 0.02857 et "less" = 0.01832 donc 20/40 offre DivGen < 8/13

Het.3.250.14 <- reintro.3ans %>% filter(scen.name %in% c("3.250.14","3.250.14-Copy","3.250.14-Copy2","3.250.14-Copy3"))
Het.3.250.14 <- na.omit(as.data.frame(tapply(Het.3.250.14$GeneDiv, Het.3.250.14$scen.name,mean)))
names(Het.3.250.14)[1] <- "mean"
Het.3.250.44 <- reintro.3ans %>% filter(scen.name %in% c("3.250.44","3.250.44-Copy","3.250.44-Copy2","3.250.44-Copy3"))
Het.3.250.44 <- na.omit(as.data.frame(tapply(Het.3.250.44$GeneDiv, Het.3.250.44$scen.name,mean)))
names(Het.3.250.44)[1] <- "mean"
Het.3.50.14 <- reintro.3ans %>% filter(scen.name %in% c("3.50.14","3.50.14-Copy","3.50.14-Copy2","3.50.14-Copy3")) 
Het.3.50.14 <- na.omit(as.data.frame(tapply(Het.3.50.14$GeneDiv, Het.3.50.14$scen.name,mean)))
names(Het.3.50.14)[1] <- "mean"
Het.3.50.44 <- reintro.3ans %>% filter(scen.name %in% c("3.50.44","3.50.44-Copy","3.50.44-Copy2","3.50.44-Copy3"))
Het.3.50.44 <- na.omit(as.data.frame(tapply(Het.3.50.44$GeneDiv, Het.3.50.44$scen.name,mean)))
names(Het.3.50.44)[1] <- "mean"

Het.5.250.14 <- reintro.5ans %>% filter(scen.name %in% c("5.250.14","5.250.14-Copy","5.250.14-Copy2","5.250.14-Copy3"))
Het.5.250.14 <- na.omit(as.data.frame(tapply(Het.5.250.14$GeneDiv, Het.5.250.14$scen.name,mean)))
names(Het.5.250.14)[1] <- "mean"
Het.5.250.44 <- reintro.5ans %>% filter(scen.name %in% c("5.250.44","5.250.44-Copy","5.250.44-Copy2","5.250.44-Copy3"))
Het.5.250.44 <- na.omit(as.data.frame(tapply(Het.5.250.44$GeneDiv, Het.5.250.44$scen.name,mean)))
names(Het.5.250.44)[1] <- "mean"
Het.5.50.14 <- reintro.5ans %>% filter(scen.name %in% c("5.50.14","5.50.14-Copy","5.50.14-Copy2","5.50.14-Copy3"))
Het.5.50.14 <- na.omit(as.data.frame(tapply(Het.5.50.14$GeneDiv, Het.5.50.14$scen.name,mean)))
names(Het.5.50.14)[1] <- "mean"
Het.5.50.44 <- reintro.5ans %>% filter(scen.name %in% c("5.50.44","5.50.44-Copy","5.50.44-Copy2","5.50.44-Copy3")) 
Het.5.50.44 <- na.omit(as.data.frame(tapply(Het.5.50.44$GeneDiv, Het.5.50.44$scen.name,mean)))
names(Het.5.50.44)[1] <- "mean"

# Diff 3 ans 44 ou 14% de surive = non
ks.test(Het.3.250.44$mean,Het.3.250.14$mean, alternative="two.sided") #Nextant p-value = 0.7714 pas de différence, GeneDiv = execo
ks.test(Het.3.50.44$mean,Het.3.50.14$mean, alternative="two.sided") #Nextant p-value = 1 et greater = 0.7788

# Diff 5 ans 44 ou 14% de surive = non
ks.test(Het.5.250.44$mean,Het.5.250.14$mean, alternative="two.sided") #p-value = 0.7714
ks.test(Het.5.50.44$mean,Het.5.50.14$mean, alternative="two.sided") #p-value = 0.2286

# Diff 3ans et 5 ans 250 ou 50ind reintroduits = oui plus important avec 250
ks.test(Het.3.250.44$mean,Het.3.50.44$mean, alternative="less") #p-value = 0.02857 différence 50 less than 250 (0.01832)
ks.test(Het.5.250.44$mean,Het.5.50.44$mean, alternative="less") #p-value = 0.02857 différence 50 less than 250 (0.01832)
ks.test(Het.3.250.14$mean,Het.3.50.14$mean, alternative="two.sided") #p-value = 0.02857 50 less than 250 (0.01832)
ks.test(Het.5.250.14$mean,Het.5.50.14$mean, alternative="two.sided") #p-value = 0.02857 50 less than 250 (0.01832)

# Pas les memes scenarios retenues avec Nall et Pextant. en commun 8 : S1, S2, S9, S10, S17, S18, S25 et S26
Nall <- allscenario_dd %>% filter(Nall>200)
dot_plot(data=Nall, project="Parameters-tests-RFGB.xlm", scenario='Basic dd',yrs=c(10,20),
         params="Nall", setcolour="Nall",plotpops='Naturelle',save2disk=FALSE)
# General graph en lignes pour densite dependance
lineplot.allscenario_dd <- line_plot_year(data=allscenario_dd, project="Parameters-tests-RFGB.xlm"
                                      ,scenario="Basic dd", params=c("Nall","Het"),plotpops='Naturelle', save2disk=FALSE)
lineplot.allscenario_dd <- line_plot_year(data=allscenario_dd, project="Parameters-tests-RFGB.xlm"
                                          ,scenario="S", params=c("Nall","Het"),plotpops='Reintroduite', save2disk=FALSE)
line_plot_year(data=best.dd.nat, project="Parameters-tests-RFGB.xlm"
                                          ,scenario="Basic dd", params=c("Nall","Het"),plotpops='Naturelle', save2disk=FALSE)
# General graph en lignes pour non densite dependance
lineplot.allscenario_nondd <- line_plot_year(data=allscenario_nondd, project="Parameters-tests-nondd-RFGB.xlm"
                                          ,scenario="S", params=c("Nall","GeneDiv"),plotpops='Naturelle', save2disk=FALSE)
lineplot.allscenario_nondd <- line_plot_year(data=allscenario_nondd, project="Parameters-tests-nondd-RFGB.xlm"
                                          ,scenario="S", params=c("Nall","GeneDiv"),plotpops='Reintroduite', save2disk=FALSE)
# Plot des valeurs de rRec positives
rrec.plot1<-ggplot(data=best.rnat.dd, aes(x=Scenario, y=rRec))+geom_point(stat="identity",aes(color=factor(rRec>0.03)), size=4)+geom_errorbar(aes(ymin=rRec,ymax=rRec+SD), show.legend=F,size=0.3)
rrec.plot1+labs(title="Moyenne et SD du taux de croissance des meilleurs scenarios des populations naturelles - simulations de 0-10 ans", x="Nom des scénarios VORTEX", y = "Moyenne du taux de croissance (rRec)")+theme_classic() 

rrec.plot2<-ggplot(data=best.rreintro.dd, aes(x=Scenario, y=rRec))+geom_point(stat="identity",aes(color=factor(rRec>0.09)), size=4)+geom_errorbar(aes(ymin=rRec,ymax=rRec+SD), show.legend=F,size=0.3)
rrec.plot2+labs(title="Moyenne et SD du taux de croissance des meilleurs scenarios des populations réintroduites - simulations de 0-10 ans", x="Nom des scénarios VORTEX", y = "Moyenne du taux de croissance (rRec)")+theme_classic() 

rrec.plot3<-ggplot(data=best.rnat.nondd, aes(x=Scenario, y=rRec))+geom_point(stat="identity",aes(color=factor(rRec>0.45)), size=4)+geom_errorbar(aes(ymin=rRec,ymax=rRec+SD), show.legend=F,size=0.3)
rrec.plot3+labs(title="Moyenne et SD du taux de croissance des meilleurs scenarios des populations naturelles - simulations de 0-10 ans", x="Nom des scénarios VORTEX", y = "Moyenne du taux de croissance (rRec)")+theme_classic() 

rrec.plot4<-ggplot(data=best.rreintro.nondd, aes(x=Scenario, y=rRec))+geom_point(stat="identity",aes(color=factor(rRec>0.2)), size=4)+geom_errorbar(aes(ymin=rRec,ymax=rRec+SD), show.legend=F,size=0.3)
rrec.plot4+labs(title="Moyenne et SD du taux de croissance des meilleurs scenarios des populations réintroduites - simulations de 0-10 ans", x="Nom des scénarios VORTEX", y = "Moyenne du taux de croissance (rRec)")+theme_classic() 

# Plot Pextinct cumulees pour les 3 populations
Pext.dd.best<-Pext.dd[[1]] %>% filter(Pext<0.009)
ggplot(data=Pext.dd[[1]], aes(x=Scenario, y=Pext))+geom_point(stat="identity", aes(color=factor(Population)), size=4)+geom_errorbar(aes(ymin=Pext-SD,ymax=Pext+SD), show.legend=F,size=0.3)
ggplot(data=Pext.dd.best[!Pext.dd.best$Population == "Metapopulation",], aes(x=Scenario, y=Pext))+geom_point(stat="identity", aes(color=factor(Population)), size=4)+geom_errorbar(aes(ymin=Pext-SD,ymax=Pext+SD), show.legend=F,size=0.3)

Pext.dd.nonbest<-Pext.nondd[[1]] %>% filter(Pext<0.009)
ggplot(data=Pext.nondd[[1]], aes(x=Scenario, y=Pext))+geom_point(stat="identity", aes(color=factor(Population)), size=4)+geom_errorbar(aes(ymin=Pext-SD,ymax=Pext+SD), show.legend=F,size=0.3)

# comparer les scenar cata avec S17/S9/S1/S25 --» pas de diff significatives avec cette variable
# extraction de allscenario_dd
selected.dd <- allscenario_dd %>% filter(scen.name %in% c("Basic dd", "S17","S9","S1","S25")) %>% filter(pop.name %in% c("Naturelle"))
cata.25 <- allcata %>% filter(scen.name %in% c("25(1)","25(2)","25(3)","25(4)")) %>% filter(pop.name %in% c("Naturelle"))
# graph de comparaison cata/non cata pop naturelle
## manque légende en bonne et du forme des scenarios (basic, et 4 replicats) et echelle y pour catastrophes
ggbase1<-ggplot(data=selected.dd, aes(x=Year, y=Nextant))+geom_ribbon(stat="identity", aes(ymin=Nextant-SD.Nextant.,ymax=Nextant+SD.Nextant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes(x=Year, y=Nextant, color=factor(scen.name)),size=1.5,show.legend=F)+labs(x="Années de simulation", y = "Taille des populations restantes naturelles")+theme_classic() 

ggcata1<-ggplot(data=cata.25, aes(x=Year, y=Nextant))+geom_ribbon(stat="identity", aes(ymin=Nextant-SD.Nextant.,
                ymax=Nextant+SD.Nextant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes
                (x=Year,y=Nextant, color=factor(scen.name)),size=1.5, show.legend=F)+theme_classic()
# Pextant 
ggplot(data=selected.dd, aes(x=Year, y=Pextant))+geom_ribbon(stat="identity", aes(ymin=Pextant-SD.PExtant.,ymax=Pextant+SD.PExtant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes(x=Year, y=Pextant, color=factor(scen.name)),size=1.5,show.legend=F)+labs(x="Années de simulation", y = "Probabilité des populations restantes")+theme_classic() 
ggplot(data=cata.25, aes(x=Year, y=PExtant))+geom_ribbon(stat="identity", aes(ymin=PExtant-SD.PExtant.,ymax=PExtant+SD.PExtant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes(x=Year, y=PExtant, color=factor(scen.name)),size=1.5,show.legend=F)+labs(x="Années de simulation", y = "Probabilité des populations restantes")+theme_classic() 

library(gridExtra)
library(cowplot)
plot_grid(ggbase1, ggcata1, labels=c("Sans assèchements", "Avec assèchements"), ncol = 2, nrow = 1)
#text(20,0.5,expression("50% threshold"),pos=3,col="plum",cex=0.95)
# scale_fill_discrete()
+labs(title="Taille des populations naturelles des meilleurs scénarios densités-dépendants sans catastrophes")

# graph de comparaison cata/non cata pop reintroduite
selected.dd2 <- allscenario_dd %>% filter(scen.name %in% c("Basic dd","S17","S25")) %>% filter(pop.name %in% c("Reintroduite"))
selected.dd3 <- allscenario_dd %>% filter(scen.name %in% c("Basicdd","S9","S1")) %>% filter(pop.name %in% c("Reintroduite"))
cata.25bis <- allcata %>% filter(scen.name %in% c("Basic-without-25","25(1)","25(4)")) %>% filter(pop.name %in% c("Reintroduite"))
cata.25biis <- allcata %>% filter(scen.name %in% c("Basic-without-25","25(2)","25(3)")) %>% filter(pop.name %in% c("Reintroduite"))

ggbase2<-ggplot(data=selected.dd2, aes(x=Year, y=Pextant))+geom_ribbon(stat="identity", aes(ymin=Pextant-SD.PExtant.,ymax=Pextant+SD.PExtant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes(x=Year, y=Pextant, color=factor(scen.name)),size=1.5,show.legend=F)+labs(x="Années de simulation", y = "Taille des populations restantes")+theme_classic() 
cata.25bis$scen.name <- fct_relevel(cata.25bis$scen.name, c("Basic-without-25","25(1)","25(4)"))
ggcata2<-ggplot(data=cata.25bis, aes(x=Year, y=PExtant))+geom_ribbon(stat="identity", aes(ymin=PExtant-SD.PExtant.,ymax=PExtant+SD.PExtant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes(x=Year, y=PExtant, color=factor(scen.name)),size=1.5,show.legend=F)+labs(x="Années de simulation", y = "Taille des populations restantes")+theme_classic() 

plot_grid(ggbase2, ggcata2, labels=c("Sans assèchements", "Avec assèchements"), ncol = 2, nrow = 1)


ggbase3<-ggplot(data=selected.dd3, aes(x=Year, y=Pextant))+geom_ribbon(stat="identity", aes(ymin=Pextant-SD.PExtant.,ymax=Pextant+SD.PExtant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes(x=Year, y=Pextant, color=factor(scen.name)),size=1.5,show.legend=F)+labs(x="Années de simulation", y = "Taille des populations restantes")+theme_classic() 
cata.25biis$scen.name <- fct_relevel(cata.25biis$scen.name, c("Basic-without-25","25(2)","25(3)"))
ggcata3<-ggplot(data=cata.25biis, aes(x=Year, y=PExtant))+geom_ribbon(stat="identity", aes(ymin=PExtant-SD.PExtant.,ymax=PExtant+SD.PExtant., color=factor(scen.name)), fill="grey", outline.type="both", show.legend=F,size=1)+geom_line(aes(x=Year, y=PExtant, color=factor(scen.name)),size=1.5,show.legend=F)+labs(x="Années de simulation", y = "Taille des populations restantes")+theme_classic() 

plot_grid(ggbase3, ggcata3, labels=c("Sans assèchements", "Avec assèchements"), ncol = 2, nrow = 1)

### Nombre d'individus reintroduits 500 vs 100
### Probabilite survie 0.44 vs 0.14
#### Assechement catastrophe avec meilleurs scenarios ####
cata25.all <- bind_rows(list(basic.cata25,cata25))
cata25.all <- cata25.all %>% filter(pop.name %in% c("Naturelle","Reintroduite"))
names(cata25.all)[20] <- "SD.GeneDiv."
dot_plot(cata25.all,project='catastrophy.xlm',scenario='Basic-25ans',yrs=c(10,25),params=c("PExtant","GeneDiv"), setcolour = "scen.name", plotpops="Naturelle", save2disk=F)
cata50.all <- bind_rows(list(basic.cata50,cata50))
cata50.all <- cata50.all %>% filter(pop.name %in% c("Naturelle","Reintroduite"))
names(cata50.all)[20] <- "SD.GeneDiv."
dot_plot(cata50.all,project='catastrophy.xlm',scenario='Basic-50ans',yrs=c(10,20,30,50),params=c("PExtant","GeneDiv"),
              setcolour = "Nextant", plotpops="Reintroduite", save2disk=F)

basic50.nat <- basic.cata50 %>% filter(basic.cata50$pop.name %in% "Naturelle")
basic50.reintro <- basic.cata50 %>% filter(basic.cata50$pop.name %in% "Reintroduite")
as.data.frame(tapply(basic50.nat$GeneDiv, basic50.nat$scen.name,mean))
as.data.frame(tapply(basic50.nat$SD.GD., basic50.nat$scen.name,mean))


cata50.nat <- cata50.all %>% filter(cata50.all$pop.name %in% "Naturelle")
cata50.reintro <- cata50.all %>% filter(cata50.all$pop.name %in% "Reintroduite")
as.data.frame(tapply(cata50.nat$GeneDiv, cata50.nat$scen.name,mean))
as.data.frame(tapply(cata50.reintro$SD.GD., cata50.reintro$scen.name,mean))

names()[1] <- "mean"

#### Assechement catastrophe avec pires scenarios ####
worst50.all <- bind_rows(list(basic.worst50,cata.worst50))
worst50.all <- worst50.all %>% filter(pop.name %in% c("Naturelle","Reintroduite"))
names(worst50.all)[20] <- "SD.GeneDiv."
dot_plot(worst50.all,project='catastrophy-worst.xlm',scenario='Basic-50ans',yrs=c(25,50),params=c("PExtant","GeneDiv"), setcolour = "scen.name", plotpops="Naturelle", save2disk=F)

worst50.nat <- worst50.all %>% filter(worst50.all$pop.name %in% "Naturelle")
worst50.reintro <- worst50.all %>% filter(worst50.all$pop.name %in% "Reintroduite")
as.data.frame(tapply(worst50.nat$GeneDiv, worst50.nat$scen.name,mean))
as.data.frame(tapply(worst50.reintro$SD.GeneDiv., worst50.reintro$scen.name,mean))

#### DATA analysis 1 rRec ####
# Calculate the mean recovery rate (Pacioni et al 2017) and compare scenarios *pop naturelle uniquement*
recov.dd.nat.10 <- rRec(allscenario_dd, project="Parameters-tests-RFGB.xlm",scenario='Basic dd', 
              ST=FALSE, runs=1000, yr0=1, yrt=10,save2disk=FALSE,dir_out='DataAnalysis/rRec')
recov.dd.nat.25 <- rRec(allscenario_dd, project="Parameters-tests-RFGB.xlm",scenario='Basic dd', 
                    ST=FALSE, runs=1000, yr0=10, yrt=25,save2disk=FALSE,dir_out='DataAnalysis/rRec')
best.rnat.dd<-recov.dd.nat.10 %>% filter(rRec>0)
recov.dd.nat.25 %>% filter(rRec>0)
best.rreintro.dd<-recov.dd.reintro.10 %>% filter(rRec>0)
recov.dd.reintro.25 %>% filter(rRec>0)

best.rnat.nondd<-recov.nondd.nat.10 %>% filter(rRec>0)
recov.nondd.nat.25 %>% filter(rRec>0)
best.rreintro.nondd<-recov.nondd.reintro.10 %>% filter(rRec>0)
recov.nondd.reintro.25 %>% filter(rRec>0)
recov.nondd.nat.10 <- rRec(allscenario_nondd, project="Parameters-tests-RFGB.xlm",scenario='Basic nondd', 
                    ST=FALSE, runs=1000, yr0=1, yrt=10,save2disk=FALSE,dir_out='DataAnalysis/rRec')
recov.nondd.nat.25 <- rRec(allscenario_nondd, project="Parameters-tests-RFGB.xlm",scenario='Basic nondd', 
                       ST=FALSE, runs=1000, yr0=10, yrt=25,save2disk=FALSE,dir_out='DataAnalysis/rRec')

# Pour les 2 populations *ne garder que reintroduite*
merge.dd.filter <- merge.dd %>% filter(pop.name==c("Reintroduite"))
merge.nondd.filter <- merge.nondd %>% filter(pop.name==c("Reintroduite"))
recov.dd.reintro.10 <- rRec(merge.dd.filter, project="Parameters-tests-RFGB.xlm",scenario='Basic dd', 
                    ST=FALSE, runs=1000, yr0=1, yrt=10,save2disk=FALSE,dir_out='DataAnalysis/rRec')
recov.dd.reintro.25 <- rRec(merge.dd.filter, project="Parameters-tests-RFGB.xlm",scenario='Basic dd', 
                    ST=FALSE, runs=1000, yr0=10, yrt=25,save2disk=FALSE,dir_out='DataAnalysis/rRec')

recov.nondd.reintro.10 <- rRec(merge.nondd.filter, project="Parameters-tests-RFGB.xlm",scenario='Basic nondd', 
                       ST=FALSE, runs=1000, yr0=1, yrt=10,save2disk=FALSE,dir_out='DataAnalysis/rRec')
recov.nondd.reintro.25 <- rRec(merge.nondd.filter, project="Parameters-tests-RFGB.xlm",scenario='Basic nondd', 
                       ST=FALSE, runs=1000, yr0=10, yrt=25,save2disk=FALSE,dir_out='DataAnalysis/rRec')
# compare SSDM avec coeff de kendall - 
#Pop Naturelle
kendal.nat.dd1 <- recov.dd.10 %>% filter(Population=="Naturelle") %>% select(SSMD)
kendal.nat.dd2 <- recov.dd.25 %>% filter(Population=="Naturelle") %>% select(SSMD)
# verifier que l ordre des facteurs de rangs calcule par SSMD sont statistiquement concordant entre les deux sets d annees 
#Pop Reintroduite
kendal.nat.dd <- recov.dd.10 %>% filter(pop.name=="Reintroduite")
kendal.nat.dd2 <- recov.dd.25 %>% filter(pop.name=="Reintroduite")
require()
kendall.DD

#### DATA analysis 2 Pextinct et regression lineaire ####
Pext.dd <- Pextinct(run.dd[[2]], project='Parameters-tests-RFGB',
        scenario='Basicdd', ST=FALSE, save2disk=FALSE)
Pext.nondd <- Pextinct(run.nondd[[2]], project='Parameters-tests-RFGB',
                    scenario='Basicnondd', ST=FALSE, save2disk=FALSE)

Pext.dd.reintro <- allscenario_dd[allscenario_dd$scen.name == c("S1","S2","S9","S10","S17","S18","S25","S26","S33","S41","S42","S49","S57","S58","S65","S66","S73","S74","S81","S82","S89","S90","S97","S105","S106","S113","S121","S122"),]
line_plot_year(data=allscenario_dd, project="Parameters-tests-RFGB"
,scenario="S17", params=c("Nall","Nextant","r.stoch"),plotpops='Reintroduite', save2disk=FALSE)

# selctionnees 
run.selected <- run.dd[[2]] %>% filter(Scenario %in% c("Basicdd", "S17","S9","S1","S25")) %>% filter(Population %in% c("Naturelle"))
run.cata.nat <- run.cata[[2]] %>% filter(Scenario %in% c("Basic-without-25","25(1)","25(2)","25(3)","25(4)")) %>% filter(Population %in% c("Naturelle"))
Pext.dd.selected <- Pextinct(run.selected, project='Parameters-tests-RFGB',
                    scenario='Basicdd', ST=FALSE, save2disk=FALSE)
Pext.cata<- Pextinct(run.cata.nat, project="catastrophy-dd", scenario="Basic-without-25",ST=FALSE, save2disk=FALSE)

run.selected2 <- run.dd[[2]] %>% filter(Scenario %in% c("Basicdd", "S17","S9","S1","S25")) %>% filter(Population %in% c("Reintroduite"))
run.cata.reintro <- run.cata[[2]] %>% filter(Scenario %in% c("Basic-without-25","25(1)","25(2)","25(3)","25(4)")) %>% filter(Population %in% c("Reintroduite"))
Pext.dd.selected2 <- Pextinct(run.selected2, project='Parameters-tests-RFGB',
                             scenario='Basicdd', ST=FALSE, save2disk=FALSE)
Pext.cata2 <- Pextinct(run.cata.reintro, project="catastrophy-dd", scenario="Basic-without-25",ST=FALSE, save2disk=FALSE)

#### GLM ####
# Remove base scenario from .run output in long format
lrun.dd.no.base <- run.dd[[2]][!run.dd[[2]]$Scenario == 'Basicdd',]
reg <- fit_regression(data=lrun.dd.no.base,
                      census=FALSE,
                      project='Parameters-tests-RFGB', scenario='S', popn=3,
                      param='N', vs=c("MortalityT"), l=2, ncand=10,
                      save2disk=FALSE)
# ne fonctionne pas, ne trouve pas GS1...

# Clean up of residual files written by glmulti
# Note, in some OS (W) these files may be locked because in use by R and have
# to be manually after the R session has been either terminated or restarted
file.remove(c('Pacioni_et_al_ST_LHS_N.modgen.back',
              'Pacioni_et_al_ST_LHS_N.mods.back'))

#### Sensitivity Test (pas utile) ####
#(not run pas avant avoir fait ST)lkup.dd <- lookup_table(data=allscenario_dd, project="~/Aurore/VORTEX10/Parameters-tests-RFGB", scenario='Parameters-tests-RFGB_S', pop=c('Reintroduite','Naturelle'), SVs=c('PExtinct','SD.PExtinct.','stoch.r', 'SD.r.', 'Nextant','SD.Nextant.','X'),save2disk=FALSE)
# Funtion collate_yr before Nadults : fichier .yr genere lorsque ST effectuÃ©e
#NadultsAllndd <- Nadults(data=allscenario_dd, gen= 3,yr0 = 1, yrt = 25, fname = "Nadults",dir_out = "DataAnalysis")

# SSMD matrix /WM : Error in sottra/sqrt(sumsq) : tableaux de tailles inadÃ©quates
ssmd_matrixdd <- SSMD_matrix(data=allscenario_dd,project="/home/aurore/Documents/Rainette_FGB/VORTEX10/Parameters-tests-RFGB"
                             ,scenario='Basic dd', 
                             params= c("PExtinct","stoch.r", "Nextant")
                             ,yrs=c(10,25), ST=FALSE, save2disk = FALSE)

  