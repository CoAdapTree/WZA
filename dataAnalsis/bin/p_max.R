## An R implementation of Mike's Pmax

pMaxTest <- function( p, setSize){
  numTests = length(p)
  sortedP = sort(p, decreasing = T)
  pThreshold = sortedP[1 + numTests - setSize] 
  sum( dbinom( setSize:numTests, numTests, pThreshold))
}

## The p_max function takes a vector of p_values (1 per species)

## Next I'll generate a set of random uniform variables as p-values to test the dist. is uniform

pMaxTest(c(0.1,0.5,0.3), 3)

pop1 <- runif(100)
pop2 <- runif(100)
pop3 <- runif(100)
pop4 <- runif(100)
pop5 <- runif(100)

M <- as.matrix( cbind(pop1,pop2,pop3, pop4, pop5) )
par( mfrow = c(2,2) )
hist( apply ( M, 1, function(x) pMaxTest(x, 5)) , main = "5 population comparison")
hist( apply ( M, 1, function(x) pMaxTest(x, 4)) , main = "4 population comparison")
hist( apply ( M, 1, function(x) pMaxTest(x, 3)) , main = "3 population comparison")
hist( apply ( M, 1, function(x) pMaxTest(x, 2)) , main = "2 population comparison")

## Yep, that checks out
par( mfrow = c(1,1) )

FiveWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x, 5)), "hochberg" )
FourWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x, 4)), "hochberg" )
ThreeWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x,3)), "hochberg" )
TwoWay <- stats::p.adjust( apply ( M, 1, function(x) pMaxTest(x,2)), "hochberg" )

p.adjust 
FiveWayHits <- FiveWay<0.05
FourWayHits <- FourWay<0.05
ThreeWayHits <- ThreeWay<0.05
TwoWayHits <- TwoWay<0.05

sum(FiveWayHits)
