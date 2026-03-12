##############################################################
#
#                  R O B U S T    M S A       R   C O D E
#
#                     AUTHOR:   URBANO LORENZO-SEVA
#                        URV, TARRAGONA (SPAIN)
#                         DATE: 29/06/2021
#
#                    EDITS: 
#                    NATASHA TONGE, GEORGE MASON UNIVERSITY, 05/03/2026
#
# DATA EXAMPLE: 61 ITEMS AND N = 578
#
# DATA TEXT FILE USED:    data_example.dat
#
# EXECUTE THE CODE AS:    source("RobustMSA_v2.r")
#
##############################################################

##############################################################
# UPDATE HERE THE NAME OF YOUR INPUT DATA FILE
# A TEXT FILE (.DAT) WITH NO LABELS, ROWNAMES, OR COLUMNS IS EXPECTED
filein <- "" 
#
# UPDATE HERE THE SHORT DESCRIPTION OF THE RUN
#
Note<- "A short descriptive note about the run"
# 
# UPDATE HERE THE NUMBER OF BOOTSTRAP SAMPLES (AT LEAST 500)
# 
K <- 3000
# 
# UPDATE HERE THE CONFIDENCE INTERVAL 
# 
C <- 0.95
#
# 
# UPDATE HERE THE MSA CRITERION TO ELIMINATE ITEMS 
# 
CRIT <- 0.50
#
##############################################################

####################  ROBUST  MSA CODE #######################
if (K<500) k <- 500
if ((C < .80) | (C > .99)) C <- 0.95  
if ((CRIT < .3) | (CRIT > .70)) CRIT <- .5  
X <- read.table(filein)
X<-na.omit(X)
N <- nrow(X)
m <- ncol(X)
R <- cor(X)
IR <-solve(R)
SS <-diag(diag(IR)^(-.5))
Q <- SS %*% IR %*% SS
RS <- R * R
QS <- Q * Q
RS <- RS-diag(diag(RS))
QS <- QS-diag(diag(QS))
MNUM <- colSums(RS)
MDEN <- colSums(RS+QS)
MSA <- MNUM/MDEN

MSAB <- matrix(0,K,m)
alive=1
alive2=1
for (i in 1:K){
    if (alive2==10) {
       cat("\014")
       alivestr <- switch(alive,"Computing... \U002B   ","Computing... \U00D7   ")
       print(alivestr)
       alive <- alive+1
       alive2 <- 1
       if (alive>2) alive <- 1
    } else alive2 <- alive2+1

    id <- floor(runif(N, min=1, max=(N+1)))
    Xi <- X[id,]
    Ri <- cor(Xi)
    IR <-solve(Ri)
    SS <-diag(diag(IR)^(-0.5))
    Q <- SS %*% IR %*% SS
    RS <- Ri * Ri
    QS <- Q * Q
    RS <- RS-diag(diag(RS))
    QS <- QS-diag(diag(QS))
    MNUM <- colSums(RS)
    MDEN <- colSums(RS+QS)
    MSAi <- MNUM/MDEN
    MSAB[i,]=MSAi
}
cat("\014")
print("Job done!   ")
print("            ")

MSAB_IC <- matrix(0,2,m)
dummy <- matrix(0,1,m)
for (i in 1:m){
    MSAB[,i] = as.matrix(sort(as.numeric(MSAB[,i]),na.last=TRUE, method ="quick"))  # Natasha Tonge, GMU edit. added as.numeric
    MSAB_IC[1,i] = MSAB[floor(K*((1-C)/2)),i]
    MSAB_IC[2,i] = MSAB[floor(K*((1+C)/2)),i]
    dummy[1,i] <- (MSAB_IC[1,i] < CRIT) 
}
eliminate <- colSums(t(dummy))

####################   END OF ROBUST MSA CODE  ######################

# ==============================================================
#  Robust MSA — Write formatted output to file (no packages)
#  Choose output format: "txt" or "html"
# ==============================================================

output_format <- "txt"          # <-- change to "txt" or use "html" output
output_file   <- "robust_msa_results" #<-- change filename to something descriptive


# --------------------------------------------------------------
# Helper: collect all output lines into a character vector
# --------------------------------------------------------------

build_lines <- function() {
  
  lines <- character(0)
  add   <- function(...) lines <<- c(lines, paste0(...))
  
  add("================================================================")
  add("")
  add("        R O B U S T   M S A   -   R   C O D E")
  add("")
  add("        Author  :  Urbano Lorenzo-Seva & Natasha Tonge")
  add("        Inst.   :  URV(Spain) & George Mason University (USA)")
  add("        Date    :  13/03/2026")
  add("")
  add("================================================================")
  add("")
  add(sprintf("  %-26s  %s",  "File in",           filein))
  add(sprintf("  %-26s  %d",  "Cases",              N))
  add(sprintf("  %-26s  %d",  "Variables",          m))
  add(sprintf("  %-26s  %d",  "Bootstrap samples",  K))
  add(sprintf("  %-26s  %s",  "Percentile",         C))
  
  if (eliminate > 0) {
    add(sprintf("  %-26s  %d", "Variables to eliminate", eliminate))
  } else {
    add(sprintf("  %-26s  %s", "Variables to eliminate", "None"))
  }
  add(sprintf("  %-26s  %s",  "File Notes",         Note))
  
  add("")
  add("----------------------------------------------------------------")
  add(sprintf("  %-6s  %-8s  %s", "Item", "MSA", "95% Confidence Interval"))
  add("----------------------------------------------------------------")
  
  for (i in 1:m) {
    flag <- if (dummy[1, i] == 1) "*" else " "
    add(sprintf("  %3d %s   %5.3f     %5.3f  --  %5.3f",
                i, flag, MSA[i], MSAB_IC[1, i], MSAB_IC[2, i]))
  }
  
  add("----------------------------------------------------------------")
  add("  * Variable flagged for elimination")
  add("")
  add("================================================================")
  
  lines
}


# --------------------------------------------------------------
# Write plain-text output
# --------------------------------------------------------------

write_txt <- function(lines, filepath) {
  con <- file(filepath, open = "wt")
  writeLines(lines, con)
  close(con)
  message("Output written to: ", filepath)
}


# --------------------------------------------------------------
# Write HTML output (no packages — pure base R string building)
# --------------------------------------------------------------

write_html <- function(lines, filepath) {
  
  is_rule    <- function(l) grepl("^={3,}|^-{3,}", trimws(l))
  is_title   <- function(l) grepl("R O B U S T", l)
  is_param   <- function(l) grepl("^  [A-Z][a-zA-Z ]+ {2,}", l)
  is_colhead <- function(l) grepl("Item.*MSA.*Confidence", l)
  is_flagged <- function(l) grepl("\\*", l) && grepl("[0-9]\\.[0-9]{3}", l)
  is_data    <- function(l) grepl("^ {2,}[0-9]+", l)
  
  row_html <- function(l) {
    esc <- gsub("&", "&amp;", l)
    esc <- gsub("<", "&lt;",  esc)
    esc <- gsub(">", "&gt;",  esc)
    
    cls <- "line"
    if (is_rule(l))    cls <- "line rule"
    if (is_title(l))   cls <- "line title"
    if (is_param(l))   cls <- "line param"
    if (is_colhead(l)) cls <- "line colhead"
    if (is_data(l))    cls <- "line data"
    if (is_flagged(l)) cls <- "line data flagged"
    
    sprintf('    <div class="%s">%s</div>', cls, esc)
  }
  
  body_rows <- paste(sapply(lines, row_html), collapse = "\n")
  
  html <- paste0('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Robust MSA Results</title>
  <style>
    body {
      background: #ffffff;
      color: #000000;
      font-family: "Courier New", Courier, monospace;
      font-size: 13px;
      line-height: 1.8;
      padding: 40px;
    }

    .line         { white-space: pre; display: block; color: #000000; }
    .line.title   { font-weight: bold; }
    .line.colhead { font-weight: bold; }
    .line.flagged { font-weight: bold; }
  </style>
</head>
<body>
', body_rows, '
</body>
</html>')
  
  con <- file(filepath, open = "wt")
  writeLines(html, con)
  close(con)
  message("Output written to: ", filepath)
}


# --------------------------------------------------------------
# Run
# --------------------------------------------------------------

lines <- build_lines()

if (output_format == "txt") {
  write_txt(lines, paste0(output_file, ".txt"))
} else {
  write_html(lines, paste0(output_file, ".html"))
}
