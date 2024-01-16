args <- commandArgs()

ratio <- read.table(args[6], header = TRUE);

png(filename = paste(args[6], ".chrscale_fast-seq.png", sep = ""), width = 2280, height = 218,
    units = "px", pointsize = 20, bg = "white", res = NA)
par(mar = c(4, 0, 0, 0))
count <- 1
widths <- vector(length = 24)
if ('chrY' %in% ratio$Chr) {
  chr_iter_list <- c(1:22, "X", "Y")
} else {
  chr_iter_list <- c(1:22, "X")
}
for (i in chr_iter_list) {
  ch <- which(ratio$Chr == paste("chr", i, sep = ""))
  widths[count] <- max(ratio$End[ch])
  count <- count + 1
}
# widths

tt <- which(ratio$Z_score > 20)
ratio$Z_score[tt] <- 20
tt <- which(ratio$Z_score < -20)
ratio$Z_score[tt] <- -20
nf <- layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24), 1, 24,
                    byrow = TRUE), widths = widths)
for (i in c(1:22, "X", "Y")) {
  tt <- which(ratio$Chr == paste("chr", i, sep = ""))
  if (length(tt) > 0) {
    plot(ratio$Start[tt], ratio$End[tt], type = "n", xlim = c(0, max(ratio$End[tt])), ylim = c(-20, 20),
         yaxt = "n", xaxt = "n", xlab = paste("chr", i))
    tt <- which(ratio$Chr == paste("chr", i, sep = "") &
                  ratio$Z_score < 3 &
                  ratio$Z_score > -3)
    segments(ratio$Start[tt], ratio$Z_score[tt], ratio$End[tt], ratio$Z_score[tt], col = colors()[88], lwd = 5)
    tt <- which(ratio$Chr == paste("chr", i, sep = "") & ratio$Z_score > 3)
    segments(ratio$Start[tt], ratio$Z_score[tt], ratio$End[tt], ratio$Z_score[tt], col = colors()[136], lwd = 5)
    tt <- which(ratio$Chr == paste("chr", i, sep = "") & ratio$Z_score < -3)
    segments(ratio$Start[tt], ratio$Z_score[tt], ratio$End[tt], ratio$Z_score[tt], col = colors()[461], lwd = 5)
    abline(h = 3)
    abline(h = -3)
  }

}

dev.off()
