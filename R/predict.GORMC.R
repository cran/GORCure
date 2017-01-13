predict.GORMC <-
function(object,...){
     arg<-list(...)
     M<-length(object$ParEst$Eta)
     P<-length(object$ParEst$Beta)

     if(is.null(arg$len)) arg$len<-100
     if(is.null(arg$new.z)) arg$new.z<-c(1,rep(0,M-1))
     if(is.null(arg$new.x)) arg$new.x<-rep(0,P)
     new.z<-as.vector(arg$new.z)
     new.x<-as.vector(arg$new.x)
     mdata<-object$mdata
     ti<-unique(c(0,na.omit(mdata$Li),na.omit(mdata$Ri)))
     if(is.null(arg$tp)){
	tp<-seq(0,max(ti),length.out=arg$len)
     }else{
	tp<-arg$tp
     }
     pzi<-exp(sum(object$ParEst$Eta*new.z))/(1+exp(sum(object$ParEst$Eta*new.z)))
     exb<-exp(sum(object$ParEst$Beta*new.x))
     Het.est1<-t(Ispline(tp[tp<=max(ti)],order=object$ParEst$order,knots=object$ParEst$knots))%*%object$ParEst$gl
     obj<-smooth.spline(tp[tp<=max(ti)],Het.est1)
     Het.est2<-predict(obj,x=tp[tp>max(ti)],deriv=0)$y

     if(object$ParEst$r>0)  sut<-(1+object$ParEst$r*c(Het.est1,Het.est2)*exb)^(-1/object$ParEst$r)
     if(object$ParEst$r==0)  sut<-exp(-c(Het.est1,Het.est2)*exb)

     surv<-1-pzi+pzi*sut
     pred<-list(CureRate=1-pzi,Survival=data.frame(Time=tp,SurvProb=surv))
     class(pred)<-"GORMC"
     class(pred)<-"predict.GORMC"
return(pred)
}
