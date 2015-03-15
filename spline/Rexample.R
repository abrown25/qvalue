library(qvalue)
x <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
y <- c(1, 0.445501, 0.439427, 0.432638, 0.431452, 0.419355, 0.418433, 0.415881, 0.409677, 0.403519, 0.394194, 0.397133, 0.393548, 0.385253, 0.378495, 0.370323, 0.362903, 0.350538, 0.348387)

df <- 3

w <- NULL
spar = NULL;
all.knots = FALSE;
nknots = .nknots.smspl;
keep.data = TRUE;
cv = FALSE;
df.offset = 0;
penalty = 1;
tol = 1e-06 * IQR(x)

contr.sp <- list(low = -1.5, high = 1.5, tol = 1e-04, eps = 2e-08,
		     maxit = 500, trace = getOption("verbose"))

xy <- xy.coords(x, y)
y <- xy$y
x <- xy$x

w <- rep_len(1, n)

ux <- x

ox <- TRUE
tmp <- cbind(w, w * y, w * y^2)[order(x), ]

wbar <- w
ybar <- y
yssw <- 0

xbar <- (x - x[1])/(x[n] - x[1])

nknots <- nknots(n)

knot <- c(rep(xbar[1], 3), xbar[seq.int(1, n, length.out = nknots)], rep(xbar[n], 3))
nk <- nknots + 2L
ispar <- -1L

spar <- as.double(0)

icrit <- 1L
dofoff <- df.offset
if (!missing(df)) {
if (df > 1 && df <= nx) {
if (!missing(cv))
warning("specified both 'df' and 'cv'; will disregard the latter")
crit <- 3L
dofoff <- df
} else warning("you must supply 1 < df <= n,  n = #{unique x} = ",
	     nx)
}
iparms <- c(icrit = icrit, ispar = ispar, iter = as.integer(contr.sp$maxit))
keep.stuff <- FALSE
ans.names <- c("coef", "ty", "lev", "spar", "parms", "crit",
	       "iparms", "ier", if (keep.stuff) "scratch")
fit <- .Fortran(C_rbart, as.double(penalty), as.double(dofoff),
			 x = as.double(xbar), y = as.double(ybar), w = as.double(wbar),
			 ssw = as.double(yssw), as.integer(nx), as.double(knot),
			 as.integer(nk), coef = double(nk), ty = double(nx), lev = double(if (is.na(cv)) 1L else nx),
			 crit = double(1), iparms = iparms, spar = spar, parms = unlist(contr.sp[1:4]),
			 scratch = double(17L * nk + 1L), ld4 = 4L, ldnk = 1L,
			 ier = integer(1L))[ans.names]
