
# Egon Pearson test function ## https://rpubs.com/seriousstats/epcs_test
# script by Tom Baguley https://seriousstats.wordpress.com/2019/09/05/chi-square-and-the-egon-pearson-correction/
########################################################################################################################
epcs.test <- function(data.cells, z.adjust.method=c('none','hommel')){
  ucs.test <- suppressWarnings(stats::chisq.test(data.cells, simulate.p.value=FALSE, correct = FALSE))
  N <- sum(data.cells) ; nrows <- dim(data.cells)[1] ; ncols <- dim(data.cells)[2]
  corrected.stat <- ucs.test$stat[[1]] * (N-1)/N
  pval <- pchisq(corrected.stat, ucs.test$par, lower.tail = FALSE)
  p.resids <- ucs.test$resid ; as.resids <- ucs.test$stdres
  pr.pv <- pchisq(ucs.test$resid^2,1, lower.tail=F)
  pr.apv <- matrix(as.matrix(p.adjust(pchisq(p.resids^2, 1, lower.tail=F), z.adjust.method[1])), nrows, ncols, byrow=F)
  asr.apv <- matrix(as.matrix(p.adjust(pchisq(as.resids^2,1, lower.tail=F), z.adjust.method[2])), nrows, ncols, byrow=F)
  output.list <- list(
    'Uncorrected Pearson Chi-Square'= ucs.test$stat,
    'Egon Pearson Chi-Square' = corrected.stat,
    'df'=prod(dim(data.cells)-1), 'p'=pval,
    'Smallest expected value (should be greater than 1)' = min(ucs.test$expected),
    'Raw data' = data.cells,
    'Pearson (standardized) residuals (z)' = p.resids,
    'two-sided p (for z)' = pr.pv, 'adjusted p (for z)' = pr.apv,
    'Chi-square contribution per cell' = p.resids^2,
    'Adjusted standardized residuals' = as.resids,
    'adjusted p (for ASRs)' = asr.apv,
    'Pearson / ASR p adjustment' =z.adjust.method)
  return(output.list)
}


########################################################################################################################
# Tests Eagon Pearson Chi2 tests pour comparer la distribution des genomes avec TET vs sans TET en catégorie d'hôte pour les catégories gut
# Sans les unknown
########################################################################################################################
 

contingence_tableH <- read.table("contingence_table_host.txt", sep="\t")


# Comparison all categories of tet(W) + IME_Rho_tet32 vs genomes without tetW/tet32
epcs.test(contingence_tableH[,c(1,2,3,4,7)])
#X-squared = 173.0635, df = 20, p-value < 3e-26, smallest expected value = 2.59

#$`Adjusted standardized residuals`
#               Ime_RhoTetW   TETW.TIR TETW.alone Ime_RhoTet32 Genomes.no.Tet
#human            3.7202009 -2.3475921 -0.1161144     6.218635      -2.247913
#cattle          -1.1859655 -0.2723463 -1.6998541    -1.638922       2.211867
#chicken          2.8775860  5.9725327 -1.2625562    -1.755650      -4.049139
#pig              0.6724246  1.9316196  5.1662708    -2.606179      -3.455129
#other domestic  -3.7688227  2.0284214 -0.6053031    -3.110979       1.778599
#other           -4.1010825 -3.9313994 -1.5992751    -2.646214       6.312535
#
#$`adjusted p (for ASRs)`
#            [,1]         [,2]         [,3]         [,4]         [,5]
#[1,] 0.004379415 2.195190e-01 9.075619e-01 1.454344e-08 2.458177e-01
#[2,] 0.907561857 9.075619e-01 5.498171e-01 6.073770e-01 2.697587e-01
#[3,] 0.072131523 6.540754e-08 9.075619e-01 5.498171e-01 1.233752e-03
#[4,] 0.907561857 4.272518e-01 6.447888e-06 1.438713e-01 1.155063e-02
#[5,] 0.003608425 3.401380e-01 9.075619e-01 3.729367e-02 5.271388e-01
#[6,] 0.001028055 1.942416e-03 6.585570e-01 1.302374e-01 7.960558e-09


# -> IME_RhoTetW : human p=0.0043, other_domestic p=0.003 et other p=0.001
# -> TetW+TIR : chicken p<10-7, other p=0.0019
# -> TetW seul : pig p<10-5
# -> IME_RhoTet32 : human p<10-7, other domestic p=0.037

