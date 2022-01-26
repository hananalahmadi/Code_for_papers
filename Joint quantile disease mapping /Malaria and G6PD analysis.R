nb <- poly2nb(map_2)
nb2INLA("map.adj", nb)
g <- inla.read.graph(filename = "map.adj")

# The separate quantile model for malaria

alpha = 0.2
r.m <- inla(formula = map_2$M.cases ~ 1+ offset(log(map_2$M.E))+
              f(u,model = "bym2", scale.model = T,
                graph = g)   ,
            family = "poisson",
            # offset = log(map_2$M.E),
            #E = map_2$M.E,
            data =  merg.data,
            control.predictor = list(compute = T),
            control.compute = list(dic = T, waic = T, cpo = T),
            control.family = list(control.link = list(model = "quantile",
                                                      quantile = alpha)),
            verbose = TRUE)




# The separate quantile model for g6pd

alpha = 0.8
r.d <- inla(formula = map_2$D.cases ~ 1+ offset(log(map_2$D.E))+
              f(u,model = "bym2", scale.model = T,
                graph = g)   ,
            family = "poisson",
            data =  merg.data,
            control.predictor = list(compute = T),
            control.compute = list(dic = T, waic = T, cpo = T),
            control.family = list(control.link = list(model = "quantile",
                                                      quantile = alpha)),
            verbose = TRUE)







# The joint quantile model

b = length(merg.data$M.cases)
y1 = c(merg.data$M.cases , rep(NA,b))
y2 = c(rep(NA,b) , merg.data$D.cases)
E1 = c(merg.data$M.E,rep(NA,b))
E2 = c(rep(NA,b) , merg.data$D.E)

yy = cbind(y1,y2)
EE = cbind(E1, E2)
mu = c(rep(1, b) , rep(2,b))
b1 = c(1:b,rep(NA ,b))

b2= c(rep(NA,b),1:b)

besagproper = c(1:b , rep(NA ,b))

Besag.c = c(rep(NA , b), 1:b)

d = data.frame(yy, mu , besagproper , Besag.c,b1,b2)  
m = as.factor(mu)

formula = yy ~ -1 +m + offset(log(E1))+ offset(log(E2))+ f(besagproper, model = "besagproper",graph=g)+
  f(Besag.c, copy="besagproper", hyper = list(beta = list(fixed = FALSE)))+
  f(b1 , model = "bym2", graph=g, scale.model=TRUE)+
  f(b2, model = "bym2", graph=g, scale.model=TRUE)




alpha = 0.2
r2 <- inla(formula,
           family = c("poisson","poisson"),
           
           control.family = list(list(control.link = list(model = "quantile",quantile = alpha)),
                                 
                                 list(control.link = list(model = "quantile",
                                                          quantile = 1-alpha))), 
           # E = EE,
           data = d,
           verbose = F,
           control.compute = list(dic = T, waic = T, cpo = T),
           control.predictor = list(compute = T))



