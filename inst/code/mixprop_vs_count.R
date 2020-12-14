count_vs_loading_plot <- function (fit, counts, s, bins = 10, title = NULL) {
  k    <- ncol(fit$L)
  bins <- seq(0,1,length.out = bins)
  bins[1] <- -1
  pdat <- NULL
  for (i in 1:k) {
    y <- tapply(counts/s,cut(fit$L[,i],bins),mean)
    pdat <- rbind(pdat,
                  data.frame(x = bins[-1],
                             y = y,
                             k = i))
  }
  pdat <- transform(pdat,k = factor(k))
  pdat <- subset(pdat,!is.na(y))
  return(ggplot(pdat,aes(x = x,y = y,color = k)) +
         geom_line(size = 0.75) +
         scale_color_manual(values = topic_colors) +
         labs(x = "mixture proportion",y = "UMI count/total",color = "topic",
              title = title) +
         theme_cowplot(font_size = 10) +
         theme(plot.title = element_text(size = 10,face = "plain")))
}

s  <- rowSums(counts)
p1 <- count_vs_loading_plot(fit_multinom,counts[,"ENSG00000105369"],s,
                            title = "CD79A")
p2 <- count_vs_loading_plot(fit_multinom,counts[,"ENSG00000105374"],s,
                            title = "NKG7")
p3 <- count_vs_loading_plot(fit_multinom,counts[,"ENSG00000184451"],s,
                            title = "CCR10")
print(plot_grid(p1,p2,p3,nrow = 1))
