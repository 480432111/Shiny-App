library(BiocManager)
library(affy)
library(GEOquery)
library(tidyverse)
library(limma)
library(shiny)
load(file = "kidneytop100genes.RData")
load(file ="rejection_labels.RData")
ui = fluidPage(
    
    sliderInput(inputId = "Cv_Folds", label = "How many CV Folds would you like?",value = 10, min = 1, max = 10),
    sliderInput(inputId = "Repeat_Value",label = "How many repeats would you like?",   value= 25, min = 1, max = 25), 
    sliderInput(inputId = "K_Value", label = "Please select a value for K for KNN", value=5, min = 1, max =10),
    plotOutput("boxplot")
    
)

server = function(input, output){
    output$boxplot <-renderPlot({
        largevar = apply(kidneytop100genes, 1, var)
        ind = which(largevar > quantile(largevar, 0.9))
        
        X = as.matrix(t(kidneytop100genes[ind,]))
        y = rejection_labels
        
        cvK = input$Cv_Folds # number of CV folds
        cv_50acc5_knn = cv_50acc5_svm = cv_50acc5_rf = c()
        cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
        
        n_sim = input$Repeat_Value ## number of repeats
        for (i in 1:n_sim) {
            
            cvSets = cvTools::cvFolds(nrow(X), cvK)  # permute all the data, into 5 folds
            cv_acc_knn = cv_acc_svm = cv_acc_rf = c()
            
            for (j in 1:cvK) {
                test_id = cvSets$subsets[cvSets$which == j]
                X_test = X[test_id, ]
                X_train = X[-test_id, ]
                y_test = y[test_id]
                y_train = y[-test_id]
                
                ## KNN
                fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = input$K_Value)
                cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
                
                ## SVM
                svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
                fit <- predict(svm_res, X_test)
                cv_acc_svm[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
                
                ## RandomForest
                rf_res <- randomForest::randomForest(x = X_train, y = as.factor(y_train))
                fit <- predict(rf_res, X_test)
                cv_acc_rf[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))
            }
            cv_50acc5_knn <- append(cv_50acc5_knn, mean(cv_acc_knn))
            cv_50acc5_svm <- append(cv_50acc5_svm, mean(cv_acc_svm))
            cv_50acc5_rf <- append(cv_50acc5_rf, mean(cv_acc_rf))
        } ## end for
        
        
        
        boxplot(list(SVM = cv_50acc5_svm, KNN = cv_50acc5_knn , RF= cv_50acc5_rf), main = "Top 100 genes - CV Accuracy for KNN, SVM and RandomForrest")
        
    })
}
shinyApp(server = server, ui = ui)