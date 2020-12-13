n <- nrow(fit$L)
k <- ncol(fit$L)

# CD79A
pdat <- data.frame(x = c(fit_multinom$L),
                   y = rep(counts[,"ENSG00000105369"],k),
                   k = factor(rep(1:k,each = n)))
p1 <- ggplot(pdat,aes(x = x,y = y,color = k)) +
  geom_smooth(method = "loess",n = 40,se = FALSE,size = 0.75) +
  scale_color_manual(values = topic_colors) +
    labs(x = "mixture proportion",y = "UMI count",
         color = "topic",title = "CD79A") +
    theme_cowplot(font_size = 10) +
    theme(plot.title = element_text(size = 10,face = "plain"))

# GATA3
pdat <- data.frame(x = c(fit_multinom$L),
                   y = rep(counts[,"ENSG00000105374"],k),
                   k = factor(rep(1:k,each = n)))
p2 <- ggplot(pdat,aes(x = x,y = y,color = k)) +
  geom_smooth(method = "loess",n = 40,se = FALSE,size = 0.75) +
  scale_color_manual(values = topic_colors) +
    labs(x = "mixture proportion",y = "UMI count",
         color = "topic",title = "NKG7") +
    theme_cowplot(font_size = 10) +
    theme(plot.title = element_text(size = 10,face = "plain"))

# GATA3
pdat <- data.frame(x = c(fit_multinom$L),
                   y = rep(counts[,"ENSG00000107485"],k),
                   k = factor(rep(1:k,each = n)))
p3 <- ggplot(pdat,aes(x = x,y = y,color = k)) +
  geom_smooth(method = "loess",n = 40,se = FALSE,size = 0.75) +
  scale_color_manual(values = topic_colors) +
    labs(x = "mixture proportion",y = "UMI count",
         color = "topic",title = "GATA3") +
    theme_cowplot(font_size = 10) +
    theme(plot.title = element_text(size = 10,face = "plain"))

print(plot_grid(p1,p2,p3,nrow = 1))
