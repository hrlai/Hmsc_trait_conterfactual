`Hmsc` does not have a function to make counterfactual predictions for new species with new traits. Here's a very crude code for it, but not tested so please proceed with caution.

# Example usage

```
library(Hmsc)

# source code
source("code/linpred_new_trait.R")

# run an example model
hM <- Hmsc(
    Y = TD$Y,
    XData = TD$X,
    XFormula = ~ x1,
    TrData = TD$Tr,
    TrFormula = ~ T1, 
    phyloTree = TD$phy
)
hM <- sampleMcmc(hM, samples = 1000, nChains = 1)

# new counterfactual trait and environmental gradients
Tr_new <- cbind(Intercept = 1,
                T1 = c(-1, 1))
X_new <- cbind(Intercept = 1,
               x1 = seq(-2, 2, 0.1))

# predict
out <- linpred_new_trait(hM, Tr = Tr_new, X = X_new)

# posterior mean of Beta of new species with new traits
Beta_mean <- apply(simplify2array(out$Beta), 1:2, mean)
Beta_mean

# posterior mean of linear predictor of new species with new traits
L_mean <- apply(simplify2array(out$L), 1:2, mean)
matplot(X_new[, -1], L_mean)

```
