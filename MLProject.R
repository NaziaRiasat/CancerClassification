#install libraries
library(readr)
library(plyr)
library(dplyr)
library(janitor)
library(float)
library(FactoMineR)
library(tidyverse)
library(stringr)
library(heatmaply)
library(factoextra)
library(ggfortify)
library(dummies)
library(caret)
library(caTools)
library(e1071)
library(entropy)
library(naivebayes)
library(MASS)
library(klaR)
library(ggpubr)


##Explore the files

##reading train data file
cancer_train <- read_csv("data_set_ALL_AML_train.csv")
#View(cancer_train)
#typeof(cancer_train)

cancer_train$call<- revalue(cancer_train$call,c("A"= 1))
cancer_train$call<- revalue(cancer_train$call,c("P"= 2))
cancer_train$call<- revalue(cancer_train$call,c("M"= 3))
head(cancer_train$call)

cancer_train1 = subset(cancer_train, select = -c(call,call_1,call_2,call_3,call_4,
                                                 call_5,call_6,call_7,call_8,call_9,call_10,
                                                 call_11,call_12,call_13,call_14,call_15,
                                                 call_16,call_17,call_18,call_19,call_20,
                                                 call_21,call_22,call_23,call_24,call_25,
                                                 call_26,call_27,call_28,call_29,call_30,
                                                 call_31,call_32,call_33,call_34,call_35,
                                                 call_36,call_37))

NewNames <- c('Gene Description', 'Gene Accession Number', '1', '2', '3', '4', '5','6',
                                                          '7','8', '9', '10','11', '12',
                                                          '13', '14', '15', '16', '17',
                                                          '18', '19', '20', '21', '22',
                                                          '23', '24', '25','26', '27',
                                                          '28', '29', '30', '31', '32', 
                                                          '33', '34', '35', '36', '37', '38')

names(cancer_train1) <- NewNames

cancer_train1.t<- t(cancer_train1)

cancer_train1.t1 <- cancer_train1.t[-c(1),] 
#typeof(cancer_train1.t1)
#summary(cancer_train1.t1)
#View(cancer_train1.t1)

# Transform first row into label
colnames(cancer_train1.t1) <- as.character.numeric_version(cancer_train1.t1[1, ], header=TRUE)
cancer_train1.t1 <- cancer_train1.t1[-1, ]

#View(cancer_train1.t1)
#typeof(cancer_train1.t1)
#class(cancer_train1.t1)
#dim(cancer_train1.t1)

#Conversion of character values into double
cancer_train1.t1 <- apply(cancer_train1.t1, 2, as.numeric)
#View(cancer_train1.t1)
#typeof(cancer_train1.t1)


#Calculating the mean & variance of each individual
mean_counts <- apply(cancer_train1.t1, 1, mean)
variance_counts <- apply(cancer_train1.t1, 1, var)

#Apply same scaling to data values
scaled_train =  scale(cancer_train1.t1 - mean_counts) / sqrt(variance_counts)
#View(scaled_train)
#typeof(scaled_train)
#str(scaled_train)

#reading independent data file
cancer_test <- read_csv("data_set_ALL_AML_independent.csv")
#head(cancer_test)
#summary(cancer_test)

cancer_test$call<- revalue(cancer_test$call,c("A"= 1))
cancer_test$call<- revalue(cancer_test$call,c("P"= 2))
cancer_test$call<- revalue(cancer_test$call,c("M"= 3))
head(cancer_test$call)

cancer_test1 = subset(cancer_test,    select = -c(call,call_1,call_2,call_3,call_4,
                                                 call_5,call_6,call_7,call_8,call_9,call_10,
                                                 call_11,call_12,call_13,call_14,call_15,
                                                 call_16,call_17,call_18,call_19,call_20,
                                                 call_21,call_22,call_23,call_24,call_25,
                                                 call_26,call_27,call_28,call_29,call_30,
                                                 call_31,call_32,call_33))

#change the column names
NewNames1 <- c('Gene Description', 'Gene Accession Number','39', '40', '41', '42', '43', '44', '45', '46',
       '47', '48', '49', '50', '51', '52', '53',  '54', '55', '56', '57', '58', '59',
       '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72')

names(cancer_test1) <- NewNames1

#take transpose
cancer_test1.t<- t(cancer_test1)
#View(cancer_test1.t)

#Remove the first row
cancer_test1.t1 <- cancer_test1.t[-c(1),] 

dim(cancer_test1.t1)
#View(cancer_test1.t1)


#To change the first row as label
colnames(cancer_test1.t1) <- as.character(cancer_test1.t1[1, ], unique=TRUE)
cancer_test1.t1 <- cancer_test1.t1[-1, ]
#View(cancer_test1.t1)
#typeof(cancer_test1.t1)

#Conversion of character values into double
cancer_test1.t1 <- apply(cancer_test1.t1, 2, as.numeric)
#View(cancer_test1.t1)

#dim(cancer_test1.t1)
#describe(cancer_test1.t1)
#typeof(cancer_test1.t1)


#Calculating the mean & variance of each individual
mean_counts1 <- apply(cancer_test1.t1, 1, mean)
variance_counts1 <- apply(cancer_test1.t1, 1, var)


#Apply same scaling to data values
scaled_test =  scale(cancer_test1.t1 - mean_counts1) / sqrt(variance_counts1)
#View(scaled_test)
#dim(scaled_test)
#typeof(scaled.test)
#str(scaled_test)


#merge both data sets
mydata=rbind(cancer_train1.t1,cancer_test1.t1)
#View(mydata)
#dim(mydata)
#head(mydata)
#tail(mydata)
#typeof(mydata)

#Histogram of un scaled data 
hist(mydata, main="Histogram of the Combined data",
     xlab="UnScaled Data",border="blue", 
     xlim=c(-10000,30000),
     col="darkmagenta",
     freq=FALSE)
lines(density(mydata))


#merge scaled data sets
mydata1=rbind(scaled_test,scaled_train)
#View(mydata1)
#dim(mydata1)
#head(mydata1)
#tail(mydata1)

#Histogram of normalize(scaled) data
hist(mydata1, main="Histogram of the scaled combined data",
     xlab="Scaled Data",border="blue", 
     col="darkmagenta",
     freq=FALSE, prob=TRUE, las=1)
lines(density(mydata1))

#removing duplicating columns
duplicated.columns <- duplicated(t(mydata1))
new.matrix <- mydata1[, !duplicated.columns]
#View(new.matrix)
#dim(new.matrix)

# Split the data into train & test
n = nrow(mydata1)
trainIndex = sample(1:n, size = round(0.75*n), replace=FALSE)
train = mydata1[trainIndex ,]
test = mydata1[-trainIndex ,]

#View(train)
#View(test)
#dim(train)
#dim(test)

#Principal Component

#Apply principal component analysis on train data
data_train.pca<- prcomp(train)# it scale the variable
summary(data_train.pca)#good results

data_train.pca1<- PCA(train,scale.unit = TRUE, ncp = 5, graph = TRUE)
summary(data_train.pca1) # better results with scaled data


#apply principal component analysis on test data
data_test.pca <- prcomp(test)# it scale the variable
summary(data_test.pca)


data_test.pca1<- PCA(test,scale.unit = TRUE, ncp = 5, graph = TRUE)
summary(data_train.pca1) # better results with scaled data

## Graphical representation of train data

#compute standard deviation of each principal component
std_dev.train <- data_train.pca$sdev

#compute variance
pr_var.train <- std_dev.train^2

#check variance of first 10 components
pr_var.train[1:10]

#proportion of variance explained
prop_varex.train <- pr_var.train/sum(pr_var.train)
prop_varex.train[1:20]

#Screeplot of principle components
plot(prop_varex.train, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")#best fit of scree plot

#cumulative scree plot
plot(cumsum(prop_varex.train), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")#good graph


fviz_eig(data_train.pca1)#best fit of scree plot


fviz_pca_biplot(data_train.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

autoplot(data_train.pca, data = train ,color=" ", label = TRUE, label.size = 3, loadings = TRUE, shape = FALSE,loadings.colour = 'blue')
# working

#3D graph

scores.train = as.data.frame(data_train.pca$x) 


library(rgl)
plot3d(scores.train[,1:3], col=c(1:4), size=10, type='p', 
       xlim = c(-20000,20000), ylim=c(-20000,20000), zlim=c(-1000,1000))
text3d(scores.train[,1]+2, scores.train[,2]+10, scores.train[,3]+2,
       texts=c(rownames(scores.train)), cex= 0.7, pos=3)

# Contributions of variables to PC1(train data)
fviz_contrib(data_train.pca, choice = "var", axes = 1, top = 10)

# Contributions of variables to PC2
fviz_contrib(data_train.pca, choice = "var", axes = 2, top = 10)


##Graphical representation for test data:


#compute standard deviation of each principal component
std_dev.test <- data_test.pca$sdev

#compute variance
pr_var.test <- std_dev.test^2

#check variance of first 10 components
pr_var.test[1:10]

#proportion of variance explained
prop_varex.test <- pr_var.test/sum(pr_var.test)
prop_varex.test[1:20]

#Screeplot of principle components
plot(prop_varex.test, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")#best fit of scree plot

#cumulative scree plot
plot(cumsum(prop_varex.test), xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b")#good graph


fviz_eig(data_test.pca)#best fit of scree plot


autoplot(data_test.pca, data =test  ,color=" ", label = TRUE, label.size = 3, loadings = TRUE, shape = FALSE,loadings.colour = 'blue')
#working

#3D graph

scores.test = as.data.frame(data_test.pca$x) 
library(rgl)
plot3d(scores.test[,1:3], col=c(1:4), size=10, type='p', 
       xlim = c(-20000,20000), ylim=c(-20000,20000), zlim=c(-1000,1000))
text3d(scores.test[,1]+2, scores.test[,2]+10, scores.test[,3]+2,
       texts=c(rownames(scores.test)), cex= 0.7, pos=3)


# Contributions of variables to PC1
fviz_contrib(data_test.pca, choice = "var", axes = 1, top = 10)

# Contributions of variables to PC2
fviz_contrib(data_test.pca, choice = "var", axes = 2, top = 10)



# Naive Bayes theorem

#Combined Actual & mydata1(merged scaled data file)

#reading actual data file
actual <- read_csv("actual.csv")

#dim(mydata1)
#View(mydata1)

mydata.12=cbind(actual,mydata1)
#View(mydata.12)# correct form of data

mydata.21 <- as.data.frame((mydata.12))
#View(mydata.21)


#Split data file into train & test with cancer column
n = nrow(mydata.21)
trainIndex = sample(1:n, size = round(0.75*n), replace=FALSE)
train11 = mydata.21[trainIndex ,]
test11 = mydata.21[-trainIndex ,]

#View(train11)
#View(test11)
#dim(train11)
#dim(test11)
#typeof(train11)
#typeof(test11)
#attributes(train11)



#transform matrix into data.21 frame
#Train data
train.df1 <- as.data.frame((train11))
#View(train.df1)
#rownames(train.df1)
#dim(train.df1)
#typeof(train.df1)

# change type of patient to character
#class(train.df1$patient)	
train.df1$patient <- as.character(train.df1$patient )

# change type of cancer to character
#class(train.df1$cancer)	
train.df1$cancer <- as.character(train.df1$cancer )


#Test data
test.df1 <- as.data.frame((test11))
#View(test.df1)
#rownames(test.df1)
#colnames(test.df1)
#dim(test.df1)
#typeof(test.df1)

# change type of patient to character
class(test.df1$patient)	
test.df1$patient <- as.character(test.df1$patient )

# change type of cancer to character
class(test.df1$cancer)	
test.df1$cancer <- as.character(test.df1$cancer )



#Naive Bayes Classification to whole data set(Best results)

newNBclassifier=naive_bayes(cancer ~ .,usekernel=FALSE, usepoisson = FALSE, mydata.21,na.action = na.pass)
newNBclassifier
summary(newNBclassifier)


# Visualize class conditional distributions corresponding to the first predictor
# with customized settings
plot(newNBclassifier, which = 1, ask = TRUE, prob = "conditional",
     arg.num = list(col = 1:3, lty = 2,
                    main = "Naive Bayes Plot", legend.position = "topleft",
                    legend.cex = 0.8))

# Visualize class marginal distributions corresponding to the first predictor
# with customized settings
plot(newNBclassifier, which = 1, ask = FALSE, prob = "marginal",
     arg.num = list(col = 1:3, lty = 1,
                    main = "Naive Bayes Plot", legend.position = "topright",
                    legend.cex = 0.55))

# Visualize class marginal distribution corresponding to the predictor "new"
# with custom colours
plot(newNBclassifier, which = 1, arg.cat = list(color = gray.colors(3)))


modelPred <- predict(newNBclassifier, mydata.21)

cMatrix <- table(modelPred, mydata.21$cancer)
cMatrix
plot(cMatrix)

confusionMatrix(cMatrix)




#Naive Bayes Classification to training data set

newNBclassifier.train=naive_bayes(cancer ~ .,usekernel=FALSE, usepoisson = FALSE, train.df1,na.action = na.pass)
newNBclassifier.train
summary(newNBclassifier.train)

# Visualize class conditional distributions corresponding to the first predictor
# with customized settings
plot(newNBclassifier.train, which = 1, ask = TRUE, prob = "conditional",
     arg.num = list(col = 1:3, lty = 2,
     main = "Naive Bayes Plot", legend.position = "topleft",
     legend.cex = 0.8))

# Visualize class marginal distributions corresponding to the first predictor
# with customized settings
plot(newNBclassifier.train, which = 1, ask = FALSE, prob = "marginal",
     arg.num = list(col = 1:3, lty = 1,
     main = "Naive Bayes Plot", legend.position = "topleft",
     legend.cex = 0.55))

# Visualize class marginal distribution corresponding to the predictor "new"
# with custom colours
plot(newNBclassifier.train, which = 1, arg.cat = list(color = gray.colors(3)))

modelPred.train <- predict(newNBclassifier.train, train.df1)

cMatrix.train <- table(modelPred.train, train.df1$cancer)
cMatrix.train
plot(cMatrix.train)

confusionMatrix(cMatrix.train)


#Naive Bayes Classification to test data set

newNBclassifier.test=naive_bayes(cancer ~ .,usekernel=FALSE, usepoisson = FALSE, test.df1, na.action = na.pass)
newNBclassifier.test
summary(newNBclassifier.test)


# Visualize class conditional distributions corresponding to the first predictor
# with customized settings
plot(newNBclassifier.test, which = 1, ask = TRUE, prob = "conditional",
     arg.num = list(col = 1:3, lty = 2,
                    main = "Naive Bayes Plot", legend.position = "topleft",
                    legend.cex = 0.8))

# Visualize class marginal distributions corresponding to the first predictor
# with customized settings
plot(newNBclassifier.test, which = 1, ask = FALSE, prob = "marginal",
     arg.num = list(col = 1:3, lty = 1,
                    main = "Naive Bayes Plot", legend.position = "topright",
                    legend.cex = 0.55))

# Visualize class marginal distribution corresponding to the predictor "new"
# with custom colours
plot(newNBclassifier.test, which = 1, arg.cat = list(color = gray.colors(3)))

modelPred.test <- predict(newNBclassifier.test, test.df1)

cMatrix.test <- table(modelPred.test, test.df1$cancer)
cMatrix.test
plot(cMatrix.test)

confusionMatrix(cMatrix.test)


#K-means clustering Algorithm

#full data

# Compute k-means with k = 5
set.seed(123)
res.km <- kmeans(mydata1, 5)
# K-means clusters showing the group of each individuals
res.km$cluster

#Graphical representation

fviz_cluster(res.km, mydata1,
             palette = c("#2E9FDF", "#00AFBB", "#E8B800" , "#3E9FDF", "#00BFBB"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#train data

# Compute k-means with k = 4
set.seed(123)
res.km.train <- kmeans(train, 4, nstart = 25)
# K-means clusters showing the group of each individuals
res.km.train$cluster

#graphical representation
fviz_cluster(res.km.train, train,
             palette = c("#2E9FDF", "#00AFBB", "#E8B800" , "#3E9FDF"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)


#test data

# Compute k-means with k = 3
set.seed(123)
res.km.test <- kmeans(test, 3, nstart = 25)
# K-means clusters showing the group of each individuals
res.km.test$cluster

#graphical representation
fviz_cluster(res.km.test, test,
             palette = c("#2E9FDF", "#00AFBB", "#E8B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)


------------
title: ML-Term-Project
Author: Nazia_Riasat 
Date: 05/05/2021 
Output: Pdf_document
df_print: kable
fontsize: 12
width: 300
matrix(runif(100), ncol = 20)
----------
        
           
options(width = 300)
matrix(runif(100), ncol = 20)


















