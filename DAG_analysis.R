library(ggdag)
library(dagitty)
library(lavaan)
library(CondIndTests)
library(dplyr)
library(GGally)
library(tidyr)
library(MKdescr)


#implied conditional independences
dataset <- read.csv("D:/data07-23.csv")
dataset <- dataset[,7:28]
str(dataset)

dataset <- dataset %>%
  select(-c(17, 19, 20, 21))

#complete cases
dataset <- dataset[complete.cases(dataset), ] 
dataset_zscore <- dataset

# zscore columns 4 to 29
for (i in 3:18) {
  dataset_zscore[[i]] <- zscore(dataset[[i]], na.rm = TRUE)
}

#DAG 
dag <- dagitty('dag {
Year -> SST12
Year -> SST3
Year -> SST34
Year -> SST4
Year -> SOI
Year -> Equatorial_SOI
Year -> NATL
Year -> SATL
Year -> TROP
Year -> CPOLR
Year -> Winds_Equator


Month -> SST12
Month -> SST3
Month -> SST34
Month -> SST4
Month -> SOI
Month -> Equatorial_SOI
Month -> NATL
Month -> SATL
Month -> TROP
Month -> CPOLR
Month -> Winds_Equator



SST12 -> SST3
SST12 -> SST34
SST12 -> SST4
SST12 -> SOI
SST12 -> Equatorial_SOI
SST12 -> NATL
SST12 -> SATL
SST12 -> TROP
SST12 -> CPOLR
SST12 -> Winds_Equator


SST3 -> SST34
SST3 -> SST4
SST3 -> SOI
SST3 -> Equatorial_SOI
SST3 -> NATL
SST3 -> SATL
SST3 -> TROP
SST3 -> CPOLR
SST3 -> Winds_Equator


SST34 -> SST4
SST34 -> SOI
SST34 -> Equatorial_SOI
SST34 -> NATL
SST34 -> SATL
SST34 -> TROP
SST34 -> CPOLR
SST34 -> Winds_Equator


SST4 -> SOI
SST4 -> Equatorial_SOI
SST4 -> NATL
SST4 -> SATL
SST4 -> TROP
SST4 -> CPOLR
SST4 -> Winds_Equator


SOI -> Equatorial_SOI
SOI -> NATL
SOI -> SATL
SOI -> TROP
SOI -> CPOLR
SOI -> Winds_Equator


Equatorial_SOI -> NATL
Equatorial_SOI -> SATL
Equatorial_SOI -> TROP
Equatorial_SOI -> CPOLR
Equatorial_SOI -> Winds_Equator


NATL -> SATL
NATL -> TROP
NATL -> CPOLR
NATL -> Winds_Equator


SATL -> TROP
SATL -> CPOLR
SATL -> Winds_Equator


TROP -> CPOLR
TROP -> Winds_Equator


CPOLR -> Winds_Equator


SST12 -> Temperature
SST3 -> Temperature
SST34 -> Temperature
SST4 -> Temperature
SOI -> Temperature
Equatorial_SOI -> Temperature
NATL -> Temperature
SATL -> Temperature
TROP -> Temperature
CPOLR -> Temperature
Winds_Equator -> Temperature


Year -> Temperature
Year -> Rain
Month -> Temperature
Month -> Rain

SST12 -> Temperature
SST3  -> Temperature
SST34 -> Temperature
SST4 -> Temperature
SOI -> Temperature
Equatorial_SOI -> Temperature
NATL -> Temperature
SATL -> Temperature
TROP -> Temperature
CPOLR -> Temperature
Winds_Equator -> Temperature


SST12 -> Rain
SST3  -> Rain
SST34 -> Rain
SST4 -> Rain
SOI -> Rain
Equatorial_SOI -> Rain
NATL -> Rain
SATL -> Rain
TROP -> Rain
CPOLR -> Rain
Winds_Equator -> Rain


SST12 -> vectors
SST3  -> vectors
SST34 -> vectors
SST4 -> vectors
SOI -> vectors
Equatorial_SOI -> vectors
NATL -> vectors
SATL -> vectors
TROP -> vectors
CPOLR -> vectors
Winds_Equator -> vectors


Temperature -> vectors 
Rain -> vectors 

Temperature -> MPI
Rain -> MPI
vectors -> MPI

Year -> excess
Month -> excess

SST12 -> excess
SST3 -> excess
SST34 -> excess
SST4 -> excess
SOI -> excess
Equatorial_SOI -> excess
NATL -> excess
SATL -> excess
TROP -> excess
CPOLR -> excess
Winds_Equator -> excess

vectors -> excess

Temperature -> excess
Rain -> excess
Rain -> Temperature

MPI -> excess

}')  


dagitty::coordinates(dag) <-  list(x=c(excess=1.5, Temperature=-0.5, SOI=-1.6, Equatorial_SOI=-1.7, SST12=-1.8, SST3=-1.9, SST4=-2.00, SST34=-2.1, NATL=-2.2, SATL=-2.3, TROP=-2.4,
                                       CPOLR=-2.5, Winds_Equator=-2.6,
                                       Year=-3.9, Month=-3.6,
                                       Rain=-0.6, 
                                       vectors=0.4, MPI=0.8),
                              
                                  y=c(excess=0.5, Temperature=0.5, SOI=1.1, Equatorial_SOI=1.2, SST12=1.3, SST3=1.4, SST4=1.5, SST34=1.6, NATL=1.7, SATL=1.8, TROP=1.9,
                                      CPOLR=2.0, Winds_Equator=2.1,
                                      Year=-3.3, Month=-3.2,
                                      Rain=2.0,
                                      vectors=-1.8, MPI=-1.5))

plot(dag)

## check whether any correlations are perfect (i.e., collinearity)
dataset_zscore %>% drop_na()

myCov <- cov(dataset_zscore)
round(myCov, 2)

myCor <- cov2cor(myCov)
noDiag <- myCor
diag(noDiag) <- 0
any(noDiag == 1)

## if not, check for multicollinearity (i.e., is one variable a linear combination of 2+ variables?)
det(myCov) < 0
## or
any(eigen(myCov)$values < 0)


## conditional independences
independences <- impliedConditionalIndependencies(dag)
print(independences)
corr <- lavCor(dataset_zscore)
summary(corr)

# plot 
r <- localTests(dag, sample.cov=corr, sample.nobs=nrow(dataset_zscore), max.conditioning.variables=3)
print(r)
plotLocalTestResults(r, xlim=c(-1,1))


# Identification
simple_dag <- dagify(
  excess ~  SST12 + SST3 + SST4 + SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month + Temperature + Rain + MPI,
  
  Temperature ~ SST12 + SST3 + SST4 + SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month, 
  
  SST12 ~ SST3 + SST4 + SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month, 
  
  SST3 ~ SST4 + SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month,
  
  SST4 ~ SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month,
  
  SST34 ~ SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month,
  
  SOI ~ Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month,
  
  Equatorial_SOI ~ NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month,
  
  NATL ~ SATL + TROP + CPOLR + Winds_Equator + Year + Month,
  
  SATL ~ TROP + CPOLR + Winds_Equator + Year + Month,
  
  TROP ~ CPOLR + Winds_Equator + Year + Month,
  
  CPOLR ~ Winds_Equator + Year + Month,
  
  Winds_Equator ~ Year + Month,
  
  Temperature ~ SST12 + SST3 + SST4 + SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month + Rain,
  
  Rain ~ SST12 + SST3 + SST4 + SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Year + Month,
  
  vectors ~ SST12 + SST3 + SST4 + SST34 + SOI + Equatorial_SOI + NATL + SATL + TROP + CPOLR + Winds_Equator + Temperature + Rain,
  
  MPI ~ Temperature + Rain + vectors,
  
  exposure = "Temperature",
  outcome = "excess",
  
  coords = list(
    x = c(
      excess = 1.5, Temperature = -0.5, SOI = -1.6, Equatorial_SOI = -1.7, SST12 = -1.8, SST3 = -1.9,
      SST4 = -2.00, SST34 = -2.1, NATL = -2.2, SATL = -2.3, TROP = -2.4,
      CPOLR = -2.5, Winds_Equator = -2.6, Year = -3.9, Month = -3.6,
      Rain = -0.6, vectors = 0.4, MPI = 0.8
    ),
    y = c(
      excess = 0.5, Temperature = 0.5, SOI = 1.1, Equatorial_SOI = 1.2, SST12 = 1.3, SST3 = 1.4,
      SST4 = 1.5, SST34 = 1.6, NATL = 1.7, SATL = 1.8, TROP = 1.9,
      CPOLR = 2.0, Winds_Equator = 2.1, Year = -3.3, Month = -3.2,
      Rain = 2.0, vectors = -1.8, MPI = -1.5
    )
  )
)


# Plot the DAG
ggdag(simple_dag, layout = "manual")
ggdag(simple_dag) + 
  theme_dag()

ggdag_status(simple_dag) +
  theme_dag()

# adjusting
adjustmentSets(simple_dag,  type = "minimal")

ggdag_adjustment_set(simple_dag, shadow = TRUE) +
  theme_dag()

