

lm.or.gam <- function(x, y, min.x, max.x) {
  # model decision tree start --------------------------------------------------------
  lin <- try(gam(y ~ x), silent = T)
  if (!inherits(lin, "try-error") & !is.na(coef(lin)[2])) {
    # try all models
    poly.2 <- try(gam(y ~ poly(x, 2)), silent = T)
    poly.3 <- try(gam(y ~ poly(x, 3)), silent = T)
    gam <- try(gam(y ~ s(x)), silent = T)
    
    model.ls <-
      list(
        lin = lin,
        poly.2 = poly.2,
        poly.3 = poly.3,
        gam = gam
      )
    # error handling, remove models that didn't work
    errors <-
      names(which(sapply(model.ls, inherits, 'try-error') == T))
    model.ls <- model.ls[!(names(model.ls) %in% errors)]
    
    # which has the smallest AIC?
    aic.df <-
      data.frame(models = as.vector(names(sapply(model.ls, AIC))),
                 AIC = as.vector(sapply(model.ls, AIC))) %>% arrange(AIC) %>% filter(!is.infinite(AIC))
    
    if (nrow(aic.df) != 0L) {
      # if any model goes through filtering, do...
      # if the best model is linear, keep the linear model
      if (aic.df[1,]$models == "lin") {
        #re-fit with lm
        model <- lm(y ~ x)
      } else {
        # if it's not linear, test against linearity
        best.model <-
          model.ls[names(model.ls) %in% aic.df[1,]$models][[1]]
        chisq.test <- anova(lin,
                            best.model, test = "Chisq")
        #aov.test <- anova(lin,
        #                  best.model)
        
        if (is.na(chisq.test$`Pr(>Chi)`[2]) |
            chisq.test$`Pr(>Chi)`[2] > 0.05) {
          # is.na(aov.test$`Pr(>F)`[2]) |
          #   aov.test$`Pr(>F)`[2] > 0.05
          
          # non-linearity is not justified
          model <- lm(y ~ x)
        } else {
          # non-linearity is justified
          # re-fit with lm if poly
          if (formula(best.model) == "y ~ poly(x, 2)") {
            model <- lm(y ~ poly(x, 2))
          } else if (formula(best.model) == "y ~ poly(x, 3)") {
            model <- lm(y ~ poly(x, 3))
          } else {
            # if gam
            model <- best.model
          }
        }
      }
    } else {
      # if there is no model that doesn't have a AIC value, return NA
      model <- NA
    }
  } else {
    # if even the linear model doesn't work return NA
    model <- NA
  }
  # model decision tree finish -------------------------------------------------------
  
  if (length(model) > 1L) {
    # if a real model was returned do...
    # we have our best model now, extract coefficients and get peak values
    # model peak extraction start -----------------------------------------------------
    if (grepl("s(", deparse(formula(model)), fixed = T) |
        grepl("poly(", deparse(formula(model)), fixed = T)) {
      # if GAM or poly then do ...
      
      # make a x-axis with a higher frequency of points to make rollapply work
      frq.x <- c(seq(min(x), max(x), by = 1), max(x))
      
      # predict y-values along high frequency x axis using the best model
      newd <- data.frame(x = frq.x)
      y.pred <- as.numeric(predict(model, newdata = newd))
      
      # apply peak detection function
      i.max <- localMaxima(y.pred)
      
      # error handling code, if no peaks fill final data frames with NA
      if (length(i.max) != 0L) {
        # if there are peaks
        # Find inflection points
        deriv3rd <- diff(diff(diff(y.pred)))
        d <- sign(deriv3rd)
        d <- cumsum(rle(d)$lengths)
        inf.pts <- d[seq.int(1L, length(d), 2L)]
        
        # Take inflection points before and after each peak
        df.ls <- list()
        
        for (i in 1:length(i.max)) {
          clos.inf <- closest(inf.pts, i.max[i])
          
          temp <-
            data.frame(
              no.peak = rep(length(i.max), times = length(clos.inf)),
              peak = rep(i, times = length(clos.inf)),
              peak.x = rep(frq.x[i.max[i]], times = length(clos.inf)),
              peak.y = rep(y.pred[i.max[i]], times = length(clos.inf)),
              closest = c("above", "below"),
              infp.x = frq.x[clos.inf],
              infp.y = y.pred[clos.inf],
              infp.slope = diff(y.pred)[clos.inf],
              stringsAsFactors = F
            )
          df.ls[[i]] <- temp
        }
        
        pk.df <- do.call(rbind, df.ls)
        
        #output
        frq.data <- data.frame(
          x.pred = frq.x,
          y.pred = y.pred,
          stringsAsFactors = F
        )
      }
    }
    # model peak extraction end ------------------------------------------------------
    
    # get model coefficients depending on the model used -----------------------------
    if (grepl("s(", deparse(formula(model)), fixed = T)) {
      # gather final output for GAM
      model.df <- data.frame(
        model = deparse(formula(model)),
        direction = ifelse(lin$coefficients[2] > 0, 1, 2),
        intercept = model$coefficients[1],
        lm.slope = lin$coefficients[2],
        x.start = min(x),
        x.end = max(x),
        mean.y = mean(y, na.rm = T),
        median.y = median(y, na.rm = T),
        sd.y = sd(y, na.rm = T),
        var.y = var(y, na.rm = T),
        dev.expl = summary(model)$dev.expl,
        adj.r.sq = summary(model)$r.sq,
        p.val = summary(model)$s.table[, 4]
      )
    } else {
      # if it's LM or poly get these....
      # same statistics as GAM just linear equivalents
      model.df <- data.frame(
        model = deparse(formula(model)),
        direction = ifelse(lin$coefficients[2] > 0, 1, 2),
        intercept = model$coefficients[1],
        lm.slope = lin$coefficients[2],
        x.start = min(x),
        x.end = max(x),
        mean.y = mean(y, na.rm = T),
        median.y = median(y, na.rm = T),
        sd.y = sd(y, na.rm = T),
        var.y = var(y, na.rm = T),
        dev.expl = summary(model)$r.squared,
        adj.r.sq = summary(model)$adj.r.squared,
        p.val = mean(summary(model)$coefficients[-1, 4]) # calculate the mean p-val across polys
      )
    }
  }
  
  # export list
  model.ls <- list(
    raw = data.frame(x = x, y = y),
    frq.data = if(exists("frq.data")) {
      frq.data
    } else {
      data.frame(x.pred = NA,
                 y.pred = NA)
    },
    model = model,
    model.df = if(exists("model.df")){
      model.df
    } else {
      data.frame(
        model = NA,
        direction = NA,
        intercept = NA,
        lm.slope = NA,
        x.start = NA,
        x.end = NA,
        mean.y = NA,
        median.y = NA,
        sd.y = NA,
        var.y = NA,
        dev.expl = NA,
        adj.r.sq = NA,
        p.val = NA
      )
    },
    pks.df = if(exists("pk.df")){
      pk.df
    } else {
      data.frame(
        no.peak = 0,
        peak = NA,
        peak.x = NA,
        peak.y = NA,
        closest = NA,
        infp.x = NA,
        infp.y = NA,
        infp.slope = NA,
        auc = NA
      )
    }
  )
  
  # return everything
  return(model.ls)
}