library(ggm)
library(cta)
library(rje)

`show_est` <- function(x){
  be <- x$beta
  cv <- x$covbeta
  df <- x$df
  se <- sqrt(diag(cv))
  w <- be/se
  pv <- 2 * pt(abs(w), df = df, lower.tail = FALSE)
  options(scipen = 5)
  x <- round(cbind(be, se, w, pv), 4)
  colnames(x) <- c("beta", "se", "wald", "p")
  x
}



# New faster powerset function

`hypercube` <- function(d) {
  U <- expand.grid(rep(list(0:1), d))
  colnames(U) <- paste("X", 1:d, sep = "")
  data.matrix(U)
}

`power_set` <- function(x) {
  sets <- lapply(
    seq_along(x),
    function(i) combn(x, i, simplify = FALSE)
  )
  unlist(sets, recursive = FALSE)
}

.P2 <- power_set(1:2)
.P3 <- power_set(1:3)
.P4 <- power_set(1:4)
.P5 <- power_set(1:5)
.P6 <- power_set(1:6)
.P7 <- power_set(1:7)
.P8 <- power_set(1:8)
.P9 <- power_set(1:9)
.P10 <- power_set(1:10)

`power_set2` <- function(x) {
  sets <- lapply(
    seq_along(x),
    function(i) t(arrangements::combinations(x, i))
  )
  lapply(sets, function(M) split(M, rep(1:ncol(M), each = nrow(M))))
}

`kronpow` <- function(C, p) {
  # Kronecker product of p equal matrices C.
  Cp <- 1
  for (i in 1:p) {
    Cp <- kronecker(Cp, C)
  }
  Cp
}

# Kronecker product of a list of matrices
`krons` <- function(...) {
  mats <- list(...)
  k <- 1
  for (mat in mats) {
    k <- k %x% as.matrix(mat)
  }
  k
}


`loglin2` <- function(d) {
  G <- 1
  K <- matrix(c(1, 1, 0, 1), 2, 2)
  for (i in 1:d) {
    G <- G %x% K
  }
  G[, -1]
}

`expTab` <- function(f, b, cod = "e") {
  ### Expand a contingency table.
  tab <- fullfact(b)
  p <- length(b)
  tabexp <- matrix(0, 0, p)
  n <- length(f)
  for (i in 1:n) {
    k <- rep(1, f[i])
    g <- tab[i, ]
    g <- t(as.matrix(as.numeric(g) - 1)) # To get the states 0, 1,...,J-1
    tabexp <- rbind(tabexp, g[k, ])
  }
  colnames(tabexp) <- names(tab)
  rownames(tabexp) <- seq_len(nrow(tabexp))
  if (cod == "e") {
    return((-1)^(tabexp + 1))
  } else if (cod == "b") {
    return(tabexp)
  }
}

`tab2data` <- function(df) {
  ### 
  tab <- xtabs(Freq ~ ., data = df)
  f <- as.vector(tab)
  d <- ncol(df)-1
  tab <- df[, 1:d]
  tabexp <- matrix(0, 0, d)
  n <- length(f)
  for (i in 1:n) {
    k <- rep(1, f[i])
    g <- tab[i, ]
    tabexp <- rbind(tabexp, g[k, ])
  }
  colnames(tabexp) <- names(tab)
  rownames(tabexp) <- seq_len(nrow(tabexp))
  tabexp
}





`fullfact` <- function(b, ch = TRUE) {
  foo <- function(x, ch) {
    u <- 1:x
    if (ch) u <- factor(u)
    u
  }
  U <- expand.grid(lapply(b, foo, ch))
  if (is.null(names(b))) {
    colnames(U) <- paste("X", seq_along(b), sep = "")
  } else {
    colnames(U) <- names(b)
  }
  U
}


`mlogit` <- function(p) {
  ### mlogit parameterization of multivariate binary table p (inv lex)
  ### b is the dimension of the table
  ### The parameters will be in the ordering of margins M
  ### i.e.  invlex order here.

  `margmat` <- function(bi, mar) {
    ### Defines the marginalization matrix
    if (mar == FALSE) {
      matrix(1, 1, bi)
    } else {
      diag(bi)
    }
  }

  `contrmat` <- function(bi, i) {
    ### Baseline contrast matrix
    if (i == FALSE) {
      1
    } else {
      cbind(-1, diag(bi - 1))
    }
  }

  d <- log(length(p), 2)
  V <- 1:d
  M <- power_set(V)
  C <- c()
  L <- c()
  for (mar in M) {
    K <- 1
    H <- 1
    for (i in V) {
      w <- is.element(i, mar)
      K <- contrmat(2, w) %x% K
      H <- margmat(2, w) %x% H
    }
    C <- blkdiag(C, K)
    L <- rbind(L, H)
  }
  eta <- C %*% log(L %*% p)
  eta <- as.vector(eta)
  names(eta) <- unlist(lapply(M, function(x) paste(x, collapse = "")))
  eta <- as.matrix(eta)
  eta
}

`marg` <- function(p, m) {
  p <- as.matrix(p)
  d <- log2(nrow(p))
  a <- array(p, rep(2, d))
  matrix(margin.table(a, m), ncol = 1)
}

margexp <- function(p, h) {
  U <- matrix(1, 2, 2)
  I <- diag(2)
  d <- log(length(p), 2)
  M <- 1
  v <- d:1
  h <- v[h] # invert because it is inv lex
  for (j in 1:d) {
    if (is.element(j, h)) {
      M <- M %x% I
    } else {
      M <- M %x% U
    }
  }
  M %*% p
}

`oddsratio` <- function(x) x[1] * x[4] / (x[2] * x[3])

# Check whether two sets of variables are jointly independent

`jointInd` <- function(p, a, b) {
  A <- marg(p, a) %*% t(marg(p, b))
  d <- as.vector(A) - marg(p, union(a, b))
  sum(abs(d))
}

`vertices` <- function(G) colnames(G)

`psetilex` <- function(n) {
  ### Ordered power set
  ### Inverse lexicographical order
  b <- rep.int(list(c(0, 1)), n)
  data.matrix(expand.grid(b))
}


#' Finds all disconnected sets.
#' G is the adjacency matrix of an undirected graph.
#' Returns a matrix.

`discon` <- function(G) {
  # G <- sign(G)
  d <- ncol(G)
  # id <- 1:d
  P <- psetilex(d)
  k <- 2^d
  ok <- c()
  for (i in 2:k) {
    sub <- as.logical(P[i, ])
    r <- sum(sub)
    Gs <- G[sub, sub, drop = FALSE]
    if (sum(Gs) != r * (r - 1)) {
      cc <- conComp(Gs)
      if (max(cc > 1)) {
        x <- sub + 0
        x[sub] <- cc
        ok <- rbind(ok, x)
      }
    }
  }
  dimnames(ok) <- list(seq_len(nrow(ok)), vertices(G))
  ok
}

#' Finds all disconnected sets.
#' G is the adjacency matrix of an undirected graph.
#' P is a list of subsets (typically the power set) to be checked.
#' Note: returns integers, not the names of the nodes.

`disconnected` <- function(G, P = power_set(id)) {
  # G <- sign(G)
  d <- nrow(G)
  id <- seq_len(d)
  # k <- length(P)
  ok <- c()
  for (i in seq_along(P)) {
    sub <- P[[i]]
    Gs <- G[sub, sub, drop = FALSE]
    A <- igraph::graph.adjacency(Gs, mode = "undirected")
    if (!igraph::is_connected(A)) {
      cc <- names(igraph::components(A)$membership)
      if (max(cc > 1)) {
        ok <- c(ok, list(cc))
      }
    }
  }
  ok
}


`moebius.inv.mat` <- function(n) {
  ### Finds the Moebius function. A matrix B  such that p = B q.
  ### n is the number of variables
  ### The order is InvLex

  b <- matrix(c(1, 0, -1, 1), 2, 2)
  B <- 1
  for (i in 1:n) {
    B <- b %x% B
  }
  B
}


`isub` <- function(v, P) {
  ### Index of all  subsets of a set v in the ordered power set P
  ### The set v is denoted by a binary vector. Ex: [1 1 0 1] = {1,2,4}

  if (missing(P)) P <- psetilex(length(v))
  x <- P <= matrix(v, nrow(P), ncol(P), byrow = TRUE)
  apply(x, 1, all)
}

`polyzero` <- function(a) {
  ### Finds the coefficients of the polynomial
  ### (a[1]-x)(a[2]-x) ... (a[n] - x)
  ### Coefficients ordered by ascending powers
  ### Is the inverse of polyroot up to ordering and scaling.

  ## Expand recursion formula
  n <- length(a)
  k <- c(1, rep(0, n))
  for (j in 1:n) {
    k[2:(j + 1)] <- k[2:(j + 1)] - a[j] * k[1:j]
  }
  rev(k)
}

`imlogit2` <- function(eta) {
  # JUNI - From a vector [logit(p1)  logit(p2)  logit2(p12) ] to a 2x2 table
  # juni(c(.2, .6, 1.2))

  p1 <- expit(eta[1])
  p2 <- expit(eta[2])
  psi <- exp(eta[3])

  if (psi == 1) {
    return(outer(c(1 - p2, p2), c(1 - p1, p1)))
  }

  a <- 1 + (p1 + p2) * (psi - 1)
  b <- -4 * psi * (psi - 1) * p1 * p2
  p11 <- 0.5 * (a - sqrt(a * a + b)) / (psi - 1)
  p01 <- p2 - p11
  p10 <- p1 - p11
  p00 <- 1 - p2 - p1 + p11
  c(p00, p10, p01, p11)
}

`imlogit` <- function(eta, print = FALSE) {
  ### Algorithm by Qaqish & Ivanova (2006)
  ### Only for binary variables: eta is in invlex order

  psi <- exp(eta)
  n <- log(length(eta), 2)
  q <- psi * 0
  Sn <- psetilex(n)
  B <- moebius.inv.mat(n)
  no <- apply(Sn, 1, sum)
  first <- (no == 1)
  q[1] <- 1
  q[first] <- psi[first] / (1 + psi[first])

  for (d in 2:n) {
    M <- Sn[no == d, , drop = FALSE]

    for (j in seq_len(nrow(M))) {
      A1 <- B[1:2^d, 1:2^d]

      last <- A1[, 2^d]
      N <- last > 0
      D <- last < 0
      m <- M[j, , drop = FALSE]

      v <- isub(m, Sn)
      qsel <- q[v][-2^d]
      a <- A1[N, -2^d] %*% qsel
      b <- A1[D, -2^d] %*% qsel

      L <- max(-a)
      U <- min(b)
      if (print) {
        cat("Margin ", (1:n)[m == 1], ": L =", L, ", U =", U, "\n")
      }
      if (L > U) {
        stop("Not strongly compatible. L > U in margin ",
          (1:n)[m == 1],
          call. = FALSE
        )
      }

      s <- 1 + m %*% 2^((1:n) - 1)

      x <- polyroot(psi[s] * polyzero(b) - polyzero(-a))

      sel <- abs(Im(x)) < 1e-10

      x <- Re(x[sel])

      pick <- (L < x) & (x < U)

      q[s] <- x[pick][1]
      # cat(x[pick], '\n')
    }
  }
  p <- B %*% q
  p
}

# Function similar to mat.mlogit from ggm
# modified for categorical data

`mat_mlogit_cat` <- function(b) {
  ### matrices of mlogit parameterization
  ### Assumes that the vectorized table is in inv_lex order.
  ### b is the dimension of the table
  ### The ordering is the same as `power_set`


  P <- power_set(seq_along(b))

  `margmat` <- function(bi, mar) {
    ### Defines the marginalization matrix
    if (mar == FALSE) {
      matrix(1, 1, bi)
    } else {
      diag(bi)
    }
  }

  `contrmat` <- function(bi, i) {
    ### Contrast matrix
    if (i == FALSE) {
      1
    } else {
      cbind(-1, diag(bi - 1))
    }
  }

  V <- seq_along(b)
  C <- c()
  L <- c()
  for (mar in P) {
    K <- 1
    H <- 1
    for (i in V) {
      w <- is.element(i, mar)
      K <- contrmat(b[i], w) %x% K
      H <- margmat(b[i], w) %x% H
    }
    C <- blkdiag(C, K)
    L <- rbind(L, H)
  }
  list(C = C, L = L)
}



#' list of contrast matrices for mlogit parameterization
#' Assumes that the vectorized table is in inv_lex order
#' The ordering is the same as `power_set`
#' @param P list of alla subsets for the responses
#' @param lev the levels of the responses
#'
#' @return a list of contrast matrices
#' @export
#'
#' @examples
`ctr_mlogit_list` <- function(P, lev) {
  `contrmat` <- function(bi, i) {
    ### Contrast matrix
    if (i == FALSE) {
      1
    } else {
      cbind(-1, diag(bi - 1))
    }
  }
  si <- 0
  V <- seq_along(lev)
  C <- vector(length(P), mode = "list")
  for (mar in P) {
    si <- si + 1
    K <- 1
    for (i in V) {
      w <- is.element(i, mar)
      K <- contrmat(lev[i], w) %x% K
    }
    C[[si]] <- K
  }
  C
}









`all_mat` <- function(lev, cc) {
  J <- length(cc)
  listC <- vector(J, mode = "list")
  listL <- vector(J, mode = "list")
  for (j in 1:J) {
    res <- lev[cc[[j]]]
    U <- mat_mlogit_cat(res)
    listC[[j]] <- U$C
    listL[[j]] <- U$L
  }
  #  list(C = listC, M = listL)
  Cbig <- do.call(blkdiag, listC)
  Mbig <- do.call(blkdiag, listL)
  list(C = Cbig, M = Mbig)
}


`Zmat` <- function(amat) {
  o <- c()
  nam <- colnames(amat)
  P <- power_set(nam) # must be sorted
  for (mar in P) {
    if (length(mar) > 1) {
      z <- conComp(amat[mar, mar])
      u <- (length(unique(z)) == 1) + 0 #       <------- correggere
    } else {
      z <- 1
      u <- 1
    }
    cat(mar, u, "\n")
    o <- c(o, u)
    #     print(paste(paste(i, collapse = ","), paste(z, collapse = " "),
    #                 paste(u==1, collapse = " "),
    #                 collapse = "/"))
  }
  Id <- diag(length(o))
  Z <- Id[, o == 1]
  Z
}

`Zmatrix` <- function(amat, b) { # b = named vector of levels
  nam <- colnames(amat) # nam = names(b) also
  P <- power_set(nam) #  sorted
  eta <- c()
  sel <- c()
  for (mar in P) {
    lab <- paste(mar, collapse = ".") # labels
    eta <- c(eta, rep(lab, prod(b[mar] - 1)))
    if (length(mar) > 1) {
      z <- conComp(amat[mar, mar])
      u <- (length(unique(z)) == 1) + 0
      new <- rep(u, prod(b[mar] - 1))
    } else {
      new <- rep(1, prod(b[mar] - 1))
    }
    sel <- c(sel, new)
  }
  names(sel) <- eta
  sel
  Id <- diag(sel)
  dimnames(Id) <- list(eta, eta)
  Z <- Id[, sel == 1]
  Z
}


`full_model` <- function(resp, regr, m, lev) {
  amat <- makeRCG(resp, regr)
  cc <- blocks(resp)
  J <- length(cc)

  len <- sum(2^unlist(lapply(cc, length)) - 1) # n. of models
  model <- rep("", len)
  #--------- Zmat -----------
  o <- c()
  jres <- c()
  for (block in cc) {
    P <- power_set(block) # Notice !
    for (i in P) {
      if (length(i) > 1) {
        z <- conComp(amat[i, i])
        u <- (length(unique(z)) == 1) + 0
      } else {
        z <- 1
        u <- 1
      }
      # print(paste(paste(i, collapse = ","), paste(z, collapse = " "), collapse = "/"))
      o <- c(o, u)
    }
    jres <- c(jres, P)
  }


  for (i in 1:len) {
    if (o[i] == 1) {
      g <- pa(jres[[i]], amat)
      if (length(g) == 0) {
        model[i] <- "~ 1"
      } else {
        model[i] <- paste("~", paste(pa(jres[[i]], amat), collapse = "*"))
      }
    } else {
      model[i] <- "~ 0"
    }
  }

  respon <- lapply(jres, paste, collapse = "*")
  block <- rep(seq_along(cc), 2^unlist(lapply(cc, length)) - 1)

  df <- data.frame(model, block, o)
  rownames(df) <- respon
  return(df[block != J, ])
  #---------------- here the model matrices are computed -------------------

  daf <- cbind(fullfact(lev), Freq = m) # the contingency table
  mat <- vector(J, mode = "list")
  for (j in 1:(J - 1)) {
    v <- unlist(cc[(j + 1):J]) # preceding blocks g_>j
    mat[[j]] <- margin_df(daf, onto = v)
  }

  # Last block is special
  mat[[J]] <- mat[[J - 1]]

  bigmat <- c()
  nam <- c()
  for (i in 1:len) {
    Z <- model.matrix(formula(df[i, 1]), data = mat[[df[i, 2]]])[drop = FALSE, ]
    print(Z[, ])
    # nl <- nlevels(mat[, resp[i]])
    # print(nl)
    ## Z <- bdiag(rep(list(Z), nl - 1))
    # bigmat <- blkdiag(bigmat, Z)
    # nam <- c(nam, colnames(Z))
  }


  colnames(bigmat) <- nam
  bigmat
}


`margin_df` <- function(daf, onto) {
  as.data.frame(margin.table(xtabs(Freq ~ ., data = daf), onto))
}

# marginalze a vector of counts m  onto the variables indexed  by j

`margin_v` <- function(m, lev, j) {
  column(margin.table(array(m, lev), j))
}

# `bdiag` <- function(li) {
#  B <- c()
#  for (M in li) {
#    B <- blkdiag(B, M)
#  }
#  B
# }

`bdiag` <- function(li) do.call(blkdiag, li)

`column` <- function(x) {
  x <- as.vector(x)
  matrix(x, length(x), 1)
}

#' makeRCG makes a regression chain graph
#'
#' @param resp
#' @param regr
#'
#' @return
#' @export
#'
#' @examples
`makeRCG` <- function(resp, regr) {
  chr <- function(x) {
    paste(as.character(x), collapse = "")
  }
  cc <- blocks(resp)
  J <- length(cc)
  mu <- lapply(resp, FUN = chr)
  if (J > 1) {
    for (j in 2:J) {
      m <- mu[[j]]
      substr(m, 1, 1) <- "+"
      mu[[j]] <- m
    }
  }
  mu <- paste(mu, collapse = "")
  bg <- UG(formula(mu))
  dg <- do.call(DAG, regr)
  G <- makeMG(bg = bg, dg = dg)
  v <- unlist(cc)
  G[v, v]
}


`blocks` <- function(lis) {
  slvars <- function(x) all.vars(x)
  lapply(lis, FUN = slvars)
}


#' makeP
#'
#' @param B a matrix of parameters
#'
#' @return P a matrix such that P %*% beta = vec(B)
#' @export
#'
#' @examples
`makeP` <- function(B) {
  vb <- as.vector(B)
  pos <- which(vb != 0)
  beta <- vb[pos]
  P <- matrix(0, prod(dim(B)), length(pos))
  P[pos, ] <- diag(length(pos))
  list(P = P, beta = beta)
}

# A che serve?

`summaryRCG` <- function(out, n, dig = 3) {
  S <- out$Shat
  I_B <- out$Bhat
  B <- diag(nrow(I_B)) - I_B
  O <- out$Ohat
  P <- makeP(B)
  beta <- P$beta
  P <- P$P
  Q <- makeP(O)
  omega <- Q$beta
  Q <- Q$P

  V <- solve(t(P) %*% (S %x% solve(O)) %*% P) / n
  se_beta <- sqrt(diag(V))

  W <- solve((1 / 2) * t(Q) %*% (solve(O) %x% solve(O)) %*% Q) / n
  se_omega <- sqrt(diag(W))
  arr <- data.frame(beta,
    se = se_beta, LB = beta - 2.58 * se_beta,
    UB = beta + 2.58 * se_beta
  )
  x <- which(B != 0, arr.ind = TRUE)
  arr_names <- paste(rownames(x), rownames(B)[x[, 2]], sep = "|")
  rownames(arr) <- arr_names
  print(round(arr, dig))
  arcs <- data.frame(omega,
    se = se_omega, LB = omega - 2.58 * se_omega,
    UB = omega + 2.58 * se_omega
  )
  x <- which(O != 0, arr.ind = TRUE)
  arc_names <- paste(rownames(x), rownames(O)[x[, 2]], sep = "~")
  rownames(arcs) <- arc_names
  print(round(arcs, dig))
  invisible()
}

#' Vec-permutation
#' Vec-permutation matrix P i.e. the matrix such that P * vec(A) = vec(A')
#' NOTE: p = ncol(A), q = nrow(A)
#' @param p
#' @param q
#'
#' @return
#' @export
#'
#' @examples
`vecper` <-
  function(p, q) {
    r <- c()
    for (i in 1:q) {
      for (j in 0:(p - 1)) {
        r <- c(r, (i + q * j))
      }
    }
    r
  }


`chainComponents` <- function(CG) {
  bg <- 0 + (CG == 100)
  con <- conComp(bg)
  v <- colnames(CG)
  d <- length(unique(con))
  cc <- vector("list", d)

  # pow <- cc
  # par <- cc

  for (j in 1:d) {
    cc[[j]] <- sort(v[con == j])
    #  pow
    #  par[[j]] <- sort(pa(cc[[j]], CG))
  }
  # list(cc = cc, parents = par)
  cc
}

`getCoords` <- function(id, d = c(0.9, 0.6)) {
  igraph::tk_fit(id)
  co <- igraph::tk_coords(id)
  ma <- apply(co, 2, max)
  shrink <- function(x, d) {
    100 * (1 - d) / 2 + d * x
  }
  co <- 100 * co %*% diag(1 / ma)
  co <- cbind(shrink(co[, 1], d[1]), shrink(co[, 2], d[2]))
  co <- round(co)
  dump("co", file = "")
}

`ones` <- function(n) {
  matrix(rep(1, n), n, 1)
}

`zeros` <- function(m, n) {
  matrix(0, m, n)
}


# The following is inefficient

`fitBGraph` <- function(amat, daf, verb = FALSE) { # same colnames
  b <- head(unlist(lapply(daf, function(x) length(unique(x)))), -1)
  flink <- function(m) {
    full_link(m,
      lev = b,
      cc = list(colnames(amat))
    )
  }
  out <- mph.fit(
    y = daf$Freq, L.fct = flink, X = Zmatrix(amat, b),
    norm.diff.conv = 1e-6,
    norm.score.conv = 1e-6, maxiter = 200, verbose = verb
  )
  cat(out$Gsq, "\n")
  cat(out$df, "d.f.\n")
  out
}


# ~~~~~~ The following function is incomplete and wrong
# ~~~~~~ Also is not used elsewhere

all_int <- function(b) { # b = named levels
  # full set of interactions eta in increasing order
  nam <- names(b)
  P <- power_set(nam)
  lp <- length(P)
  eta <- vector(lp, mode = "list")
  for (j in 1:lp) {
    mar <- P[[j]]
    inter <- paste(mar, collapse = ".")
    new <- rep(inter, prod(b[mar] - 1))
    eta[[j]] <- new
  }
  eta
}




# ~~~~~ Used by mlogit_link  ~~~ bottleneck

"contrmat" <- function(bi, i) {
  ### Contrast matrix
  if (i == FALSE) {
    1
  } else {
    cbind(-1, diag(bi - 1))
  }
}

`ctr_list` <- function(P, b) {
  g <- length(P)
  Clist <- vector(g, mode = "list")
  for (i in seq_len(g)) {
    Clist[[i]] <- ctr(P[[i]], b) # b is fixed
  }
  Clist
}




# Alternative (same efficiency)
`ctr_list2` <- function(P, b) {
  g <- length(P)
  lapply(seq_len(g), function(i) ctr(P[[i]], b))
}

# Some examples where ctr should issue an error
# > ctr(c(1,4), c(2,2,2))
# [,1] [,2]
# [1,]   -1    1
# > ctr(c(1,4), c(2,2,2,2))
# [,1] [,2] [,3] [,4]
# [1,]    1   -1   -1    1

`ctr` <- function(mar, b) {
  contrmat <- function(bi, i) {
    if (i == FALSE) {
      1
    } else {
      cbind(-1, diag(bi - 1))
    }
  }
  V <- seq_along(b) # integers not names
  C <- 1
  for (i in V) {
    w <- is.element(i, mar)
    C <- contrmat(b[i], w) %x% C
  }
  C
}

# mlogit_link_old <- function(m, b, P) {
#   tbl <- array(m, b)
#  # P <- power_set(seq_along(b))
#   g <- length(P)
#   eta <- vector(g, mode = "list")
#   for (i in seq_len(g)) {
#     C <- ctr(P[[i]], b)
#     eta[[i]] <- C %*% log(as.vector(marginTable(tbl, P[[i]])))
#   }
#   column(unlist(eta))
# }

#' Multivariate logistic link
#'
#' @param m the vectorized table of counts of all variables
#' @param lev a subset of the levels for the responses
#' @param P the power set relative to the responses
#' @param Clist the list of contrast matrices relative to the reponses.
#'
#' @return the multivariate logistic parameters for the conditional
#'         distribtioin of the responses given the rest
#' @export
#'
#' @examples
`mlogit_link` <- function(m, lev, P, Clist) {
  p <- prod(lev) # 6 if lev = c(3,2)
  k <- length(m) / p # 12/6 = 2    i
  g <- length(P) # 3 = length(Clist)            j
  eta <- vector(g, mode = "list")
  H <- sapply(0:(k - 1), function(i) i * p + 1:p)
  for (j in seq_len(g)) {
    mar <- P[[j]]
    for (i in seq_len(k)) {
      h <- H[, i]
      tabh <- array(m[h], lev)
      eta[[j]][[i]] <- Clist[[j]] %*% log(as.vector(marginTable(tabh, mar)))
    }
  }
  column(unlist(eta))
}


`linkf_cat` <- function(m, lev) {
  p <- prod(lev)
  out <- mat_mlogit_cat(lev)
  M <- out$L # fixed
  C <- out$C # fixed
  k <- length(m) / p
  eta <- matrix(0, 0, p - 1)
  for (i in 0:(k - 1)) {
    h <- i * p + 1:p
    etanew <- C %*% log(M %*% m[h])
    eta <- rbind(eta, t(etanew))
  }
  column(eta)
}

`linkf_cat2` <- function(m, lev) {
  p <- prod(lev)
  out <- mat_mlogit_cat(lev)
  M <- out$L # fixed
  C <- out$C # fixed
  k <- length(m) / p
  eta <- matrix(0, p - 1, 0)
  for (i in 0:(k - 1)) {
    h <- i * p + 1:p
    etanew <- C %*% log(M %*% m[h])
    eta <- cbind(eta, etanew)
  }
  column(eta)
}




# Variables in cc should follow the order of lev
# lev elements should be named !

`full_link` <- function(m, lev, cc) {
  J <- length(cc)
  daf <- cbind(fullfact(lev), Freq = m)
  full <- vector(J, mode = "list")
  for (j in seq_len(J)) {
    res <- lev[cc[[j]]]
    mar <- unlist(cc[j:J])
    dafm <- margin_df(daf, onto = mar)
    full[[j]] <- linkf_cat(dafm$Freq, res)
  }
  do.call(rbind, full)
}

`full_mlogit_link` <- function(m, b, cc) {
  J <- length(cc)
  tbl <- array(m, b)
  Eta <- vector(J, mode = "list")
  for (j in 1:J) {
    resp <- b[cc[[j]]]
    mar <- unlist(cc[j:J])
    mt <- marginTable(tbl, mar)
    eta <- mlogit_link(mt, resp)
    Eta[[j]] <- eta
  }
  do.call(rbind, Eta)
}

`var_lev` <- function(tbl) {
  sapply(dimnames(tbl), length)
}

# Margin table of the predecessors  SIMPLIFY

`marg_pre` <- function(tbl, cc, i) {
  k <- match(unlist(cc[i]), names(dimnames(tbl)))
  marginTable(tbl, k)
}

`array_perm` <- function(tbl, neworder) {
  k <- match(neworder, names(dimnames(tbl)))
  aperm(tbl, k)
}

`node_degree` <- function(A) {
  apply(A, 1, sum)
}

`laplacian` <- function(A) {
  diag(node_degree(A)) - A
}



#' Mlogit parametrization
#'
#' @param amat named adjacency matrix
#' @param b named vector of levels of the variables
#' @param P list of all the subsets
#'
#' @return a list
#' @export
#'
#' @examples
`par_mlogit` <- function(amat, b, P) {
  g <- length(P)
  Pnam <- lapply(P, function(v) colnames(amat)[v])
  dis <- disconnected(amat, Pnam)
  v <- sapply(ctr_list(P, b), nrow)
  w <- cumsum(v)
  par <- vector(g, mode = "list")
  for (j in seq_len(g)) {
    par[[j]] <- seq(w[j], len = v[j]) - v[j] + 1
  }
  out <- list(
    mar = P,
    inter = Pnam, par = par, id = c()
  )
  # Indices of zero-constrained interaction parameters
  id <- c()
  sel <- match(dis, out$inter)
  for (t in sel) {
    p <- out$par[t]
    id <- c(id, p[[length(p)]])
  }
  out$id <- id
  out
}


#' Label the parameters of a marginal loglinear model
#'
#' @param lis the list returned by par_mlogit

`make_labels` <- function(lis) {
  mar <- lis$mar
  inter <- lis$inter
  par <- lis$par
  lab <- c()
  for (m in seq_along(mar)) {
    A <- paste(mar[[m]], collapse = ",")
    itm <- inter[m]
    for (i in seq_along(itm)) {
      B <- paste(itm[[i]], collapse = ".")
      parmi <- par[m][[i]]
      for (j in seq_along(parmi)) {
        lab <- c(lab, paste("[", A, "] : ", B, "(", j, ")", sep = ""))
      }
    }
  }
  lab
}

#' Reorder an array
#'
#' @param tab array
#' @param ord a vector specifying the new order of the variables
#'
#' @return a permuted array
#' @export
#'
#' @examples
`areorder` <- function(tab, ord) {
  perm <- match(ord, names(dimnames(tab)))
  aperm(tab, perm)
}


`combine_factors` <- function(cc, df) {
  # df <- as.data.frame(tab)
  dat <- matrix(0, nrow(df), length(cc))
  nam <- c()
  j <- 0
  for (sub in cc) {
    j <- j + 1
    if (length(sub) > 1) {
      v <- apply(df[, sub, drop = FALSE], 1, paste,
        collapse = "."
      )
      nam <- c(nam, paste(sub, collapse = "."))
    } else {
      v <- df[, sub]
      nam <- c(nam, sub)
    }
    dat[, j] <- v
  }
  colnames(dat) <- nam
  if ("Freq" %in% colnames(df)) {
    out <- data.frame(dat, Freq = df$Freq)
  } else {
    out <- data.frame(dat)
  }
  for (j in seq_along(cc)) {
    out[, j] <- factor(out[, j])
  }
  out
}

`combine` <- function(set, data) {
  if (length(set) > 1) {
    v <- apply(data[, set, drop = FALSE], 1, paste,
      collapse = "."
    )
    v <- as.data.frame(factor(v))
    nam <- paste(set, collapse = ".")
    names(v) <- nam
    return(data.frame(nam = v, data))
  } else {
    warning("The variables should be at least 2.")
  }
}

`ciTest` <- function(a, b, c, data) {
  A <- strsplit(a, split = ".", fixed = TRUE)
  B <- strsplit(b, split = ".", fixed = TRUE)
  C <- strsplit(c, split = ".", fixed = TRUE)
  if (c == "") {
    cc <- c(A, B)
    return(ci.test(a, b, data = combine_factors(cc, df = data)))
  } else {
    cc <- c(A, B, C)
    return(ci.test(a, b, c, data = combine_factors(cc, df = data)))
  }
}


`essential_graph` <- function(dagx) {
  a <- essentialGraph(dagx)
  a[(a == 1) & (t(a) == 1)] <- 10
  a
}

`row_factors` <- function(ft) {
  ### ft: a flat table.
  ### Value: a data frame with the factors associated to the rows.
  ### Requires expand.grid from R 1.7.1
  vars <- attr(ft, "row.vars")
  k <- length(vars)
  expl <- expand.grid(vars[k:1])
  expl[, k:1, drop = FALSE]
}

`lrtest` <- function(out) {
  c(
    LRT = out$Gsq,
    df = out$df,
    p = 1 - pchisq(out$Gsq, df = out$df)
  )
}

`print_coef` <- function(out) {
  be <- round(out$beta, 3)
  be_se <- round(sqrt(diag(out$covbeta)), 2)
  pval <- round(2 * (1 - pnorm(abs(be / be_se))), 3)
  h <- cbind(be, be_se, pval)
  colnames(h) <- c("beta", "se", "p")
  h
}

# This program gives a flat table for the order D, G, J, B
# ft <- ftable(B + J ~ G + D, tab4)
# fft <- format(ft)
# fft[,1:2] <- fft[,2:1]
# noquote(fft)


`ci_bn` <- function(A, B, C=NULL) {
  if (is.null(C))
  {
    G <- diag(2)
    dimnames(G) <- list(c(A, B), c(A, B))
  }
  else{
    k <- length(C)
    G <- rbind(matrix(0, 2, 2 + k),
               cbind(matrix(1, k, 2), outer(1:k, 1:k, "<")))
    V <- c(A, B, C)
    dimnames(G) <- list(V, V)
  }
  iG <- graph_from_adjacency_matrix(G)
  as.bn(iG)
}

`ci_test` <- function(A, B, C=NULL, data) {
  V <- c(A, B, C)
  data <- data[, V]
  bn <- ci_bn(A, B, C)
  sat <- amat(bn)
  sat[A, B] <- 1
  bnsat <- as.bn(graph_from_adjacency_matrix(sat))
  
  Lhat <- score(bn,    data, type = "loglik")
  Lsat <- score(bnsat, data, type = "loglik")
  2 * (Lsat - Lhat)
}

`ci_fact` <- function (amat) 
{
  amat <- topSort(amat)
  d<- ncol(amat)
  amat <- amat[d:1,d:1]
  nod <- rownames(amat)
  ind <-vector(d-1, mode = "list")
  for (r in 1:(d-1)) {
    v <- nod[r]
    pre <- nod[-(1:r)]
    pa <-  nod[amat[, r] == 1]
    w <- setdiff(pre, pa)
    if(all(w=="")) next
    # w <- paste(w, collapse = ".")
    b <- c(v, list(w), list(pa))
    ind[[r]] <- b
  }
  ind[lapply(ind, is.null)==FALSE]
}

`ci_test2` <- function(A, B, C=NULL, data, test = "mi-g") {
  b <- length(B)
  ci <- rep(0, b)
  df <- rep(0, b)
  cin <- vector(length = b, mode = "list")
  for (i in 1:b) {
    if (is.null(C)){
      cin[[i]] <- list(A, B[i])
      te <- ci.test(A, B[i], data = data, test=test)
    }
     else{
      cin[[i]] <- list(A, B[i], C) 
      te <- ci.test(A, B[i], C, data = data, test=test)
    }
    ci[i] <- te$statistic
    df[i] <- te$parameter
    C <- union(C, B[1:i])
  }
  tests <- cbind(ci, df)
  list(stat = margin.table(tests, 2), cin = cin)
}  

`ci_AB` <- function(A, B, C=NULL) {
  b <- length(B)
  cin <- vector(length = b, mode = "list")
  for (i in 1:b) {
    if (is.null(C)){
      cin[[i]] <- list(A, B[i], C)
    }
    else{
      cin[[i]] <- list(A, B[i], C) 
    }
    C <- union(C, B[1:i])
  }
  cin
}


ci_test <- function(A,
                      B,
                      C = NULL,
                      test = test,
                      data = data) {
  b <- length(B)
  cin1 <- ci_AB(A, B, C)
  label <- NULL
  lrt <- NULL
  df <- NULL
  pval <- NULL
  
  for (j in 1:b) {
    cin2 <- ci_AB(cin1[[j]][[2]],  cin1[[j]][[1]], cin1[[j]][[3]])
    for (k in  1:length(cin2)) {
      x <- cin2[[k]][[1]]
      y <- cin2[[k]][[2]]
      z <- cin2[[k]][[3]]
      if(is.null(z)){
        indep <- paste(paste(x, "_|_"), paste(y))
        m <- ci.test(x, y, test = test, data = data)
      }
      else {
        indep <- paste(paste(x, "_|_"), paste(y, "|") ,
                       paste(z, collapse = ", "))
        m <- ci.test(x, y, z, test = test, data = data)
      }
      label <- c(label, indep)
      lrt <- c(lrt, m$statistic)
      df <- c(df, m$parameter)
      pval <- c(pval, m$p.value)
    }
  }

out <- data.frame(cond_indep = label, test = lrt, df)  
outr <- data.frame(cond_indep = label, test = round(lrt, 4), 
                  df = round(df))
marg <- c("Sum", round(apply(out[,2:3], 2, sum),4))
rbind(outr, marg)
}  

`ci_print` <- function(ci_list) {
  for (k in 1:length(ci_list)) {
    x <- ci_list[[k]][[1]]; x <- paste(x, collapse = " ")
    y <- ci_list[[k]][[2]]; y <- paste(y, collapse = " ")
    z <- ci_list[[k]][[3]]; z <- paste(z, collapse = " ")
    if (is.null(z)) {
      out <- paste(paste(x, "_|_"), paste(y))
    }
    else{
      out <- paste(paste(x, "_|_"),
                   paste(y, "|"),
                   paste(z, collapse = ", "))
    }
    cat(out, "\n")
  }
}

## TODO

`ci_test_DAG` <- function(dag, test = test, data = data){
  ci <- ci_fact(dag)
  d <- length(ci)
  for(i in 1:d){
    x <- ci[[i]][[1]]
    y <- ci[[i]][[2]] 
    z <- ci[[i]][[3]]
    if(is.null(z)){
      print(ci_test(x, y, test, data))
    }
    else{
    print(ci_test(x, y, z, test, data))
    }
  }
}



`LRT_norm` <- function(S, Shat,n){
  Khat <- solve(Shat)
  -n * log(det(S %*% Khat)) 
}

`dev_norm` <- function(a,b,c, data){
  abc <- c(a,b,c)
  data <- data[,abc]
  S <- cov(data)
  n <- nrow(data)
  Shat <- S
  f <- S[a,c, drop = FALSE] %*% solve(S[c,c]) %*% S[c,b, drop = FALSE]
  Shat[a,b] <- f
  Shat[b,a] <- f
  Khat <- solve(Shat)
  dev <- - n * log(det(S %*% Khat))
  df <- length(a) * length(b)
  list(dev = dev, df= df)
}
