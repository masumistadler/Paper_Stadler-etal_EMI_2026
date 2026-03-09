# FT-ICR-MS ----------------------------------------------------------------------------------
# from ftmsRanalysis package, adapted to use for data frames that include H/C and O/C information
# instead of their specified peak format
# Source code at: https://github.com/EMSL-Computing/ftmsRanalysis/R
# code taken from assign_class() and getVanKrevelenCategories()
# Reference:
# Bailey, V. L., Smith, A. P., Tfaily, M., Fansler, S. J., & Bond-Lamberty, B. (2017).
# Differences in soluble organic carbon chemistry in pore waters sampled from different pore size domains.
# Soil Biology and Biochemistry, 107, 133-143.

assign_class <- function(hc, oc){
  # define chemical compound classes and their ranges
  vanKrevelenCategoryLogic <- data.frame(
    HC.low = c(">= 1.5", ">= 0.8", ">= 1.5", ">= 0.8", ">= 1.5", ">= 1.5", ">= 0.8", ">= 0.2", ">= 0"),
    HC.high = c("<= 2.5", "< 2.5", "<= 2.3", "< 1.5", "<= 2.5", "<= 2.2", "< 1.5", "< 0.8", "< Inf"),
    OC.low = c("> 0", ">= 0", "> 0.3", "> 0.125", "> 0.7", "> 0.55", "> 0.65", ">= 0", ">= 0"),
    OC.high = c("<= 0.3", "<= 0.125", "<= 0.55", "<= 0.65", "<= 1.5", "<= 0.7", "<= 1.1", "<= 0.95", "< Inf"),
    NC.low = rep(">= 0", 9),
    NC.high = rep("< Inf", 9),
    PC.low = rep(">= 0", 9),
    PC.high = rep("< Inf", 9),
    NP.low = rep(">= 0", 9),
    NP.high = rep("< Inf", 9),
    O.low = rep(">= 0", 9),
    O.high = rep("< Inf", 9),
    N.low = rep(">= 0", 9),
    N.high = rep("< Inf", 9),
    P.low = rep(">= 0", 9),
    P.high = rep("< Inf", 9),
    S.low = rep(">= 0", 9),
    S.high = rep("< Inf", 9),
    mass.low = rep(">= 0", 9),
    mass.high = rep("< Inf", 9)
  )
  # remove logical signs
  bound_match <- vanKrevelenCategoryLogic
  bound_match <- data.frame(lapply(bound_match, function(x) as.numeric(as.character(gsub("> |>= |< |<= ", "", x)))))
  category <- c("Lipid","Unsat Hydrocarbon", "Protein","Lignin","Carbohydrate","Amino Sugar","Tannin","Cond Hydrocarbon","Other")
  rownames(bound_match) <- as.character(category)
  
  df <- data.frame(hc, oc) %>% data.table::setDT()
  #apply conditions
  df[hc >= bound_match$HC.low[1] & hc <= bound_match$HC.high[1] & oc > bound_match$OC.low[1] & oc <= bound_match$OC.high[1],
     cat := "Lipid"]
  df[hc >= bound_match$HC.low[2] & hc < bound_match$HC.high[2] & oc >= bound_match$OC.low[2] & oc <= bound_match$OC.high[2],
     cat := "Unsat Hydrocarbon"]
  df[hc >= bound_match$HC.low[3] & hc <= bound_match$HC.high[3] & oc > bound_match$OC.low[3] & oc <= bound_match$OC.high[3],
     cat := "Protein"]
  df[hc >= bound_match$HC.low[4] & hc < bound_match$HC.high[4] & oc > bound_match$OC.low[4] & oc <= bound_match$OC.high[4],
     cat := "Lignin"]
  df[hc >= bound_match$HC.low[5] & hc <= bound_match$HC.high[5] & oc > bound_match$OC.low[5] & oc <= bound_match$OC.high[5], 
     cat := "Carbohydrate"]
  df[hc >= bound_match$HC.low[6] & hc <= bound_match$HC.high[6] & oc > bound_match$OC.low[6] & oc <= bound_match$OC.high[6], 
     cat := "Amino Sugar"]
  df[hc >= bound_match$HC.low[7] & hc < bound_match$HC.high[7] & oc > bound_match$OC.low[7] & oc <= bound_match$OC.high[7],
     cat := "Tannin"]
  df[hc >= bound_match$HC.low[8] & hc < bound_match$HC.high[8] & oc >= bound_match$OC.low[8] & oc <= bound_match$OC.high[8],
     cat := "Cond Hydrocarbon"]
  
  df[is.na(cat), cat := "Other"]
  
  out <- as.vector(df$cat)
  return(out)
}