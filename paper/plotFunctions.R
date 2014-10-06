library(ggplot2)

results <- read.table(file="timings", header=T)
results <- cbind(results, paste(results[, 1], results[, 2], results[, 3], sep='_'))
results <- results[order(results[, 4]),]

plot.data <- results[!duplicated(results[, 6]),]
plot.data$Software <- gsub("L", "l", plot.data$Software)
plot.data$Type = factor(plot.data$Type, levels = c("Bootstrap", "Normal"), labels = c("Bootstrap", "Standard"))

pdf("Timings.pdf")
ggplot(plot.data, aes(y = CPU, x = PValues, col=Software)) + geom_line() +
          facet_wrap(~Type) +
          theme(strip.text = element_text(size=rel(1.5)), legend.text = element_text(size=rel(1.5)),
                legend.position = "top", legend.title = element_text(size=rel(1.5)),
                axis.title = element_text(size=rel(1.5))) +
          labs(y = "CPU time (seconds)", x = "Number of p values")
dev.off()

## rm(list=ls(all=T))
## results <- read.table(file="timings", header=T)
## results <- cbind(results, paste(results[, 1], results[, 2], results[, 3], sep='_'))
## results <- results[!is.na(results[, 5]), ]
## results <- results[order(results[, 5]), ]
## plot.data <- results[!duplicated(results[, 6]), ]

## plot.data <- plot.data[order(plot.data[,1]),]
## plot.data <- plot.data[order(plot.data[,3]),]
## plot.data <- plot.data[order(a[,2]),]

## pdf("mem_use.pdf")
## ggplot(data = data.frame(PValues = plot.data[1:8, 1], Memory = plot.data[9:16,5] / plot.data[1:4,5],
          ##                   Type = rep(c("Bootstrap", "Standard"), each = 4)),
          ## aes(x=PValues, y = Memory, col=Type)) +
          ## geom_line() +
          ## labs(x = "Number of p values", y = "Ratio of memory use") +
          ## theme(legend.position = "top", legend.text = element_text(size=rel(1.5)),
          ##       axis.title = element_text(size=rel(1.5)),
          ##       legend.title = element_text(size=rel(1.5)))
## dev.off()
