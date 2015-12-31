closedpMS.t <- function(X, dfreq = FALSE, h=NULL, h.control=list(), 
                        maxorder = t - 1, forced = 1:t, stopiflong = TRUE, ...) {
  
  call <- match.call()
  
  #########  Validation des arguments en entr?e et initialisation de variables  #########
  
  # Validation des arguments et obtention de la valeur de t
  valid.one(dfreq,"logical")
  Xvalid <- valid.X(X=X, dfreq=dfreq)
  X <- Xvalid$X
  t <- Xvalid$t
  
  valid.h.out <- valid.h(h=h, m=NULL, call=call,
                                    values=c("Chao","LB","Poisson","Darroch","Gamma","Normal"))
  h <- valid.h.out$h
  htype <- valid.h.out$htype
  
  if(!is.list(h.control)) stop("'h.control' must be a list")
  theta <- valid.theta(theta = h.control$theta, htype = htype)
  neg <- valid.neg(neg = h.control$neg, htype = htype)
  initsig <- valid.initsig(initsig = h.control$initsig, htype = htype)
  method <- valid.method(method = h.control$method, htype = htype)
  # initcoef ne peut pas ?tre fourni car le nombre de param?tres change selon le mod?le

  # maxorder, forced et stopiflong sont valid?s lors de l'appel ? la fonction getAllModels
  
#   if(is.null(call$maxorder)) {
#     maxorder <- t - 1
#   } else {
#     if (!(maxorder %in% 1:(t-1))) 
#       stop("'maxorder' must be an integer between 1 and 't' - 1 inclusively")
#   }
#   
#   valid.one(stopiflong, type = "logical")
#   if(stopiflong && (t==6 && maxorder >=3 || t>=7 && maxorder >=2))
#     stop("the number of models to fit is large, therefore this command should be long to run,",
#          "\nset 'stopiflong' to FALSE if you really want to run it")
  
  #########  Cr?ation de variables communes ? tous les mod?les  #########   
  
  ### Cr?ation du vecteur de variable r?ponse Y
  Y <- histfreq.t(X=X,dfreq=dfreq)
  n <- sum(na.rm=TRUE,Y)
  
  ### Cr?ation de la variable offset
  cst <- 0  
  
  
  #########  Boucle sur tous les mod?les  ######### 
  
  # Identification de tous les mod?les ? ajuster
  name <- getAllModels(t=t, maxorder=maxorder, forced=forced, stopiflong=stopiflong)
  modeles <- getFormulaFromName(name=name)
  
  # Objets pour stocker les r?sultats
  tableau <- matrix(NA_real_, nrow=length(modeles), ncol=8)
  dimnames(tableau) <- list(name, c("abundance", "stderr", "bias",
                                    "deviance", "df", "AIC", "BIC", "infoFit"))
  fit.err <- vector(mode="list", length=length(modeles))
  fit.warn <- vector(mode="list", length=length(modeles))
  names(fit.err) <- names(fit.warn) <- name
  if (htype == "Chao") {
    neg.eta <- vector(mode="list", length=length(modeles))
    names(neg.eta) <- name
  }
  
  # Boucle qui ajuste tous les mod?les
  for (j in 1:length(modeles)) {
    
    ### Cr?ation de la matrice X (qui est propre au mod?le)
    getmX.out <- getmX(typet=TRUE, t=t, t0=t, m=NULL, h=h, theta=theta, 
                                  mX=modeles[[j]])
    mX. <- getmX.out$mX. 
    nbcap <- getmX.out$nbcap
    nca <- getmX.out$nca
    
    ### Ajustement du mod?le
    fit.out <- closedp.fitone(n = n, Y = Y, mX. = mX., nbcap = nbcap, nca = nca, 
                                         cst = cst, htype = htype, neg = neg, 
                                         initsig = initsig, method = method, ...)
    
    if(!is.null(fit.out$fit.err)) fit.err[[j]] <- fit.out$fit.err
    if(!is.null(fit.out$fit.warn)) fit.warn[[j]] <- fit.out$fit.warn
    if (htype == "Chao") neg.eta[[j]] <- fit.out$neg.eta
    
    # Calculs ? faire seulement si glm a produit une sortie
    # (on laisse dans tableau et param les NA mis lors de l'initialisation
    #  en cas d'erreur lors de l'ajustement du mod?le)
    if (is.null(fit.out$fit.err)) {
      
      # On met les statistiques d'ajustement du mod?le dans tableau
      tableau[j, 3:7] <- fit.out$resultsFit
      
      # On estime N
      intercept <- if(htype == "Normal") fit.out$fit$parameters[1, 1] else coef(fit.out$fit)[1]
      stderr.intercept <- if(htype == "Normal") fit.out$fit$varcov[1, 1] else vcov(fit.out$fit)[1, 1]
      tableau[j, 1:2] <- getN(n = n, intercept = intercept, stderr.intercept = stderr.intercept)
      
      # Avertissement pour grand biais au besoin
      biasWarn <- getBiasWarn(N = tableau[j, "abundance"], bias = tableau[j, "bias"])
      if(!is.null(biasWarn)) fit.warn[[j]] <- c(fit.warn[[j]], biasWarn) 
      
    }
    
    # code pour les conditions mis dans le tableau
    tableau[j, 8] <- getInfo(err = fit.err[[j]], warn = fit.warn[[j]])
    
  }
  
  # Sortie des r?sultats
  ans <- list(n = n, t = t, results = tableau, fit.err = fit.err, fit.warn = fit.warn)
  if (htype == "Chao") ans <- c(ans, neg.eta)
  class(ans) <- "closedpMS"
  ans
}

#' @export
print.closedpMS <- function(x, ...) {
  cat("\nNumber of captured units:",x$n,"\n\n")
  
  if (!is.null(x$results)) {
    cat("Abundance estimations and model fits for the models with the smallest BIC:\n")
    tableau <- x$results[order(x$results[, "BIC"]), ]
    tabprint(tab = tableau[1:min(10,nrow(tableau)), ], 
                        digits = c(1,1,1,3,0,3,3,NA), warn = x$fit.warn, ...)
  }
  
  cat("\n")
  invisible(x)
}

#' @export
plot.closedpMS <- function(x, main="Models comparison based on BIC", omitOutliers = TRUE, ...){
  N <- x$results[, "abundance"]
  BIC <- x$results[, "BIC"]
  
  if (omitOutliers) {
    # On retire les valeurs extr?mes : < q1 - 1.5*IQR  ou  > q3 + 1.5*IQR
    Nq1 <- quantile(N,0.25)
    Nq3 <- quantile(N,0.75)
    keepN <- N > Nq1 - 1.5*(Nq3-Nq1) & N < Nq3 + 1.5*(Nq3-Nq1)
    
    BICq1 <- quantile(BIC,0.25)
    BICq3 <- quantile(BIC,0.75)
    keepBIC <- BIC > BICq1 - 1.5*(BICq3-BICq1) & BIC < BICq3 + 1.5*(BICq3-BICq1)
    
    N <- N[keepN & keepBIC]
    BIC <- BIC[keepN & keepBIC]
  }
  
  # Graphique
  plot(x = N, y = -BIC, xlab = "abundance", ylab = "-BIC", main = main, ...)
  abline(v = N[order(BIC)][1], col= "blue") 
}


# ---------------------------------------------------------------------------- #

#' @export
getAllModels <- function(t, maxorder = t - 1, forced = 1:t, stopiflong = TRUE) {
  
  # Validations et initialisations
  if (!(t %in% 2:9)) stop("'t' must be an integer between 2 and 9 inclusively")
  
  if (!(maxorder %in% 1:(t-1))) 
    stop("'maxorder' must be an integer between 1 and 't' - 1 inclusively")
  
  forced <- unique(as.character(forced))
  if (!all(unlist(strsplit(forced, "")) %in% 1:t))
    stop("the only accepted characters in 'forced' are ", 
         paste(1:(t-1),collapse=", "), " and ", t)
  
  valid.one(stopiflong, type = "logical")
  if(stopiflong && (t==6 && maxorder >=3 || t>=7 && maxorder >=2))
    stop("the number of models to enumerate is large, ",
         "therefore this command should be long to run,",
         "\nset 'stopiflong' to FALSE if you really want to run it")
  
  # Informations sur les termes forc?s dans le mod?le
  forcedOrder <- nchar(forced)
  minorder <- if(length(forced) == 0) 1 else max(forcedOrder)
  if (minorder > maxorder) 
    stop("the order of one of the forced term is larger than 'maxorder'")
  
  # Liste qui contiendra tous les mod?les possibles
  modeles <- list()
  # (un mod?le = un vecteur des termes de haut de hi?rarchie formant son nom)
  # (l'objet modeles change de taille au fil des it?rations, car il n'y a pas de
  #  formule simple pour calculer le nombre total de mod?les)
  
  #---------------------------------------------------------------------------#
  
  # Boucle sur les ordres
  # on ins?re dans les mod?les les termes un ordre ? la fois,
  # en d?butant par les termes de l'ordre le plus ?l?v? possible, soit t-1
  for (i in maxorder:1) {
    
    # Tous les termes possibles de cet ordre
    termesPossibles <- vapply(combn(1:t,i, simplify = FALSE), 
                              paste, collapse="", FUN.VALUE = "a")
    
    # Termes de cet ordre forc?s dans le mod?le
    iforced <- forced[forcedOrder == i]
    
    # ?tape A : int?gration de termes aux mod?les pr?c?dents, un mod?le ? la fois
    if (length(modeles) > 0) {
      for (j in 1:length(modeles)) {
        
        # Identification des termes d'ordre inf?rieur qui sont obligatoirement
        # pr?sents dans le mod?le (par hi?rarchie), qui ne doivent donc pas
        # appara?tre dans le nom du mod?le
        termesInfForces <- NULL
        for (k in 1:length(modeles[[j]])) {
          combin <- combn(unlist(strsplit(modeles[[j]][k], split = "")), i, simplify = FALSE)
          termesInfForces <- c(vapply(combin, paste, collapse="", FUN.VALUE = "a"),
                               termesInfForces)
        }
        iforcedNonHierar <- iforced[!(iforced %in% termesInfForces)]
        
        # Identification des termes qui pourraient appara?tre dans le nom du mod?le
        termesInfLibres <- setdiff(termesPossibles, termesInfForces)
        termesInfLibresNotForced <- setdiff(termesInfLibres, iforced)
        
        # Boucle sur le nombre de termes d'ordre inf?rieur non forc?s 
        # int?gr?s au mod?le
        if (length(termesInfLibresNotForced) > 0) {          
          for (k in length(termesInfLibresNotForced):1) {
            
            # ?num?ration de toutes les combinaisons de termes d'ordre inf?rieur
            # libres et non forc?s
            termesInfChoisis <- combn(termesInfLibresNotForced, k, simplify = FALSE)
            
            # Ajout des termes forc?s s'il y en a et s'ils ne sont pas
            # pr?sents par hi?rarchie
            termesInfChoisis <- lapply(termesInfChoisis, c, iforcedNonHierar)
            
            # Ajout de ces combinaisons aux termes d?j? dans le mod?le
            ajout <- lapply(termesInfChoisis, c, modeles[[j]])
            # Pour remettre les termes dans l'ordre voulu
            na <- k + length(iforcedNonHierar)
            ajout <- lapply(ajout, '[', c((na+1):(na+length(modeles[[j]])), 1:na))
            
            # On ajoute ces mod?les ? la liste de tous les mod?les
            modeles <- c(modeles, ajout)
          }
        }
        
        # Ajout des termes forc?s de l'ordre i sans autres termes de l'ordre i
        if (length(iforcedNonHierar) > 0) {
          modeles <- c(modeles, list(c(modeles[[j]], iforcedNonHierar)))
        }
        
        # Il faut retirer le mod?le j de la liste si :
        # - au moins un terme forc? est dans termesInfLibres 
        #   (car le mod?le ne contient pas le terme en question)
        # - 
        if(length(termesInfLibres) > length(termesInfLibresNotForced)) { 
          modeles[[j]] <- character(0)
        }
        
      }  
    }
    
    # ?tape B : ajout de mod?les sans termes d'ordre supp?rieur ? i
    # On n'ins?re que les mod?les contenant les termes forc?s de
    # l'ordre i trait? ? cette it?ration, s'il y en a.   
    if(i >= minorder) {
      
      termesPossiblesNotForced <- setdiff(termesPossibles, iforced)
      
      # Boucle sur le nombre de termes possibles non forc?s
      if (length(termesPossiblesNotForced) > 0) {
        for (j in length(termesPossiblesNotForced):1) {
          modeles <- c(modeles, 
                       lapply(combn(termesPossiblesNotForced,
                                    j, simplify = FALSE), 
                              c, iforced))
        }
      }
      
      # Ajout du mod?le contenant uniquement les termes forc?s d'ordre i,
      # s'il y a des termes forc?s de cet ordre
      if (length(iforced) > 0) modeles <- c(modeles, list(iforced))
    }
    
    # ?tape C : Retirer les ?l?ments de la liste qui sont des 
    # vecteurs vides, s'il y a lieu
    modeles <- modeles[lapply(modeles, length) > 0]
    
  } 
  
  # Pour tous les mod?les, on met les termes formant sont nom entre crochets,
  # s?par?s par des virgules
  noms <- paste0("[", vapply(modeles, paste, collapse = ",", FUN.VALUE = "a"), "]")
  # On ajoute le mod?le contenant seulement une ordonn?e ? l'origine,
  # s'il n'y a pas de termes forc?s dans le mod?le
  if(length(forced) == 0) noms <- c(noms, "[]")
  # Sortie : on retourne les noms trouv?s
  noms
}


