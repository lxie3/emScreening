#after performing step 1 and step 2 of ttscreening, return the ttscreening passing rate for identified mediators
#use lavaan package to return alpha std_alpha p_alpha beta std_beta p_beta alpha*beta std_alpha*beta
##Start the joint screening method, adjusted by ttscreening
##input are two datasets, Phenotype, DNAm
##output are matrix of selected mediators and 16 columns

EMscreening=function(Phenotype,DNAm,Continuous=TRUE,Cutoff.Joint=0.1,Iterations=100,Train.SigLevel=0.05,Test.SigLevel=0.05,Percentage=0.5){
  
  p=dim(DNAm)[1]#number of potential mediators
  n=dim(DNAm)[2]#sample size
  class(Phenotype)="numeric"
  joint.out.all=rep(NA, p)
  
  ### step 1 :screening each cg
  for (i in 1:p){
    tempdata <- cbind(X=as.vector(Phenotype[,1]),Y=as.vector(Phenotype[,2]), M=as.vector(t(DNAm[i, ])))
    colnames(tempdata)=c("X","Y","M")
    if (Continuous==TRUE) {
      joint.out.all[i]=lm_indirect_jointest(tempdata)
    }else{
      joint.out.all[i]=glm_indirect_jointest(tempdata)
    }
  }
  
  #step 2: use 0.1 as cutoff
  Final_joint=matrix(NA,nrow=p,ncol=2)
  colnames(Final_joint)<-c("pass","tt_rate")
  rownames(Final_joint)<-rownames(DNAm)
  
  for (i in 1:p){
    if (joint.out.all[i] >=Cutoff.Joint) {
      Final_joint[i]=0
    }else{
      temp <- cbind(X=as.vector(Phenotype[,1]),Y=as.vector(Phenotype[,2]), M=as.vector(t(DNAm[i, ])))
      colnames(temp)=c("X","Y","M")
      
      sig_path=rep(NA,Iterations)
      for (k in 1:Iterations) {
        set.seed(k)
        
        sample <- sample.int(n = n, size = floor(.67*n), replace = F)## Every loop generate IDs for train and test datasets
        
        ### we need to try each
        train <- temp[sample, ]
        test  <- temp[-sample, ]
        if (Continuous==TRUE) {
          Path.out.train=lm_indirect_jointest(train)
          Path.out.test=lm_indirect_jointest(test)
        }else{
          Path.out.train=glm_indirect_jointest(train)
          Path.out.test=glm_indirect_jointest(test)
        }
        
        sig_path[k] <- ifelse (Path.out.train<Train.SigLevel &Path.out.test<Test.SigLevel,1,0)
      }
      #Final_joint[i,1]=ifelse (sum(sig_path)>Percentage*Iterations,1,0) # this Percentage*Iterations is a percentage of passing rate,default value of Percentage*Iterations = 0.5, 0<q<=1,
      Final_joint[i,1]=ifelse (sum(sig_path)>Percentage*Iterations,1,0)
      Final_joint[i,2]=sum(sig_path)/Iterations
    }
  }
  
  #Final_joint[,1]==1 will give identified mediators from the joint screening method
  #will run lavaan package for each of identified mediator
  selected<-subset(Final_joint,Final_joint[,1]==1)
  No.selected<-dim(selected)[1] #number of identified mediators after ttscreening
  
  if (No.selected==0){
    print("No mediators identified from the joint screening method")
  }else{
    
    selected_ID=rownames(selected)
    
    DNAm_selected<-DNAm[selected_ID,]
    
    
    #two datasets: Phenotype, DNAm_seleted
    
    Results<-matrix(NA,nrow=No.selected,ncol=16)
    rownames(Results)<-rownames(selected)
    colnames(Results)<-c("alpha*beta est"," alpha*beta se","alpha*beta z","alpha*beta pvalue","alpha est"," alpha se","alpha z","alpha pvalue","beta est","beta se","beta z","beta pvalue","total est","total se","total z","total pvalue")
    for (i in 1:No.selected){
      tempdata <- cbind(X=as.vector(Phenotype[,1]),Y=as.vector(Phenotype[,2]), M=as.vector(t(DNAm[i, ])))
      colnames(tempdata)=c("X","Y","M")
      tempdata<-as.data.frame(tempdata)
      model <- ' # direct effect
             Y  ~ c*X
             M ~ a*X
             Y~ b*M
             # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
      if (Continuous==TRUE) {
        fit <- lavaan::sem(model, data = tempdata)
      }else{
        fit <- lavaan::sem(model, data = tempdata, ordered ="Y")
      }
      Estimate=as.matrix(lavaan::parameterEstimates(fit))
      Results[i,1:4]<-Estimate[Estimate[,c("label")]=="ab",5:8]
      Results[i,5:8]<-Estimate[Estimate[,c("label")]=="a",5:8]
      Results[i,9:12]<-Estimate[Estimate[,c("label")]=="b",5:8]
      Results[i,13:16]<-Estimate[Estimate[,c("label")]=="total",5:8]
    }
    Final_results<-cbind(selected,Results) #this is the results table we want the user to get eventually
    Final_results #18 columns
    return(Final_results)
  }
}
