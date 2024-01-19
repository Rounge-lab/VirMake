require(ggplot2)
input <- snakemake@input[[1]]
output <- snakemake@output[[1]]

# prepare data
benchmarks <- read.csv(input)
benchmarks <- benchmarks[, -2]

# create plots
dir.create(output)
for (col in colnames(benchmarks)) {
    if (col != "rule_name") {
        ggplot(benchmarks, aes(x = rule_name, y = benchmarks[, col])) +
            geom_col() +
            xlab("") +
            ylab(col) +
            coord_flip()
        ggsave(
            paste0(output, "/", col, ".png"),
            plot = last_plot(), width = 10, height = 10
        )
    }
}
