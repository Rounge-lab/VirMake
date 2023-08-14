options(warn=-1)

require("shiny")
require("shinydashboard")
require("shinythemes")
require("tidyverse")
require("ggplot2")
require("reshape2")


args = commandArgs(trailingOnly=TRUE)
output.path <- args[1]
port <- args[2]
stats.path <- paste(output.path, "statistics", sep = "/")
combined.sample.stats <- read.table(paste(stats.path, "Combined_Sample_stats.tsv", sep = "/"), header = TRUE, sep = "\t")
sample.stats.vibrant <- read.table(paste(stats.path, "Sample_stats_vibrant.tsv", sep = "/"), header = TRUE, sep = "\t")
sample.stats.virsorter2 <- read.table(paste(stats.path, "Sample_stats_virsorter2.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.AMGs <- read.table(paste(stats.path, "vOTU_AMGs.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.mapped.to.reads <- read.table(paste(stats.path, "vOTU_mapped_to_reads.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.relative.abundance <- read.table(paste(stats.path, "vOTU_Relative_Abundance.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.stats.combined <- read.table(paste(stats.path, "vOTU_stats_combined.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.stats.vibrant <- read.table(paste(stats.path, "vOTU_Stats_vibrant.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.stats.virsorter2 <- read.table(paste(stats.path, "vOTU_Stats_virsorter2.tsv", sep = "/"), header = TRUE, sep = "\t")
compared.samples <- read.table(paste(stats.path, "compared_samples_comparisonsTable.tsv",
                               sep="/"), sep="\t", header=TRUE)

UI <- fluidPage(theme = shinytheme("cerulean"),
    navbarPage(
        "VirMake results",
        tabPanel("QC",
            titlePanel("QC results"),
            sidebarPanel(
                checkboxGroupInput(
                    inputId = "sel.sample",
                    label = "Filter by sample:",
                    choices = sort(c("All", colnames(subset(vOTU.relative.abundance, select=-ID)))),
                    selected = "All",
                ),
                checkboxGroupInput(
                    inputId = "sel.quality",
                    label = "Filter by quality:",
                    choices = sort(c("All", unique(vOTU.stats.combined$checkv_quality))),
                    selected = "All",
                ),
            ),
            mainPanel(
                plotOutput("abundance.plot", height = "1000px"),
                plotOutput("checkv.quality.plot", height = "1000px")
            )
        ),
        tabPanel("Assembly",
            titlePanel("Assembly results"),
            sidebarPanel(
                checkboxGroupInput(
                    inputId = "sel.sample.assembly1",
                    label = "Filter by sample:",
                    choices = sort(c("All", unique(combined.sample.stats$Sample))),
                    selected = "All"
                ),
            ),
            mainPanel(
                plotOutput("assembly1.plot", height = "1000px"),
            )
        ),
        tabPanel("Taxonomy",
            titlePanel("Taxonomy results"),
            sidebarPanel(
                checkboxGroupInput(
                    inputId = "sel.family",
                    label = "Filter by family:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Family)[unique(vOTU.stats.combined$Family) != "n.a."])),
                    selected = "All"
                ),
                checkboxGroupInput(
                    inputId = "sel.subfamily",
                    label = "Filter by subfamily:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Subfamily)[unique(vOTU.stats.combined$Subfamily) != "n.a."])),
                    selected = "All"
                ),
                checkboxGroupInput(
                    inputId = "sel.genus",
                    label = "Filter by genus:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Genus)[unique(vOTU.stats.combined$Genus) != "n.a."])),
                    selected = "All"
                ),
            ),
            mainPanel(
                plotOutput("family.plot", height = "1000px"),
                plotOutput("subfamily.plot", height = "1000px"),
                plotOutput("genus.plot", height = "1000px"),
            )
        ),
        tabPanel("Functional analysis",
            titlePanel("Functional analysis results"),
            sidebarPanel(
                checkboxGroupInput(
                    inputId = "sel.type",
                    label = "Filter by type:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Type)[unique(vOTU.stats.combined$Type) != "n.a."])),
                    selected = "All"
                ),
                checkboxGroupInput(
                    inputId = "sel.provirus",
                    label = "Filter by provirus:",
                    choices = sort(c("All", unique(vOTU.stats.combined$provirus))),
                    selected = "All"
                ),
            ),
            mainPanel(
                plotOutput("type.plot", height = "1000px"),
                plotOutput("provirus.plot", height = "1000px"),
            )
        ),
        navbarMenu("Advanced",
            tabPanel("DRAMv",
                titlePanel("DRAMv results"),
                includeHTML(paste(output.path, "DRAMv", "distilled", "product.html", sep = "/"))
            ),
            tabPanel("inStrain compare",
                titlePanel("inStrain compare results"),
                mainPanel(
                    plotOutput("instrain.hist", height = "1000px"),
                    plotOutput("mean.cov.overlap", height = "1000px"),
                    plotOutput("mean.bases.count", height = "1000px"),
                    plotOutput("percent.genome.compared", height = "1000px"),
                    plotOutput("mean.length", height = "1000px"),
                    plotOutput("mean.consensus.SNPs", height = "1000px"),
                    plotOutput("mean.pop.SNPs", height = "1000px"),
                )
            )
        )
    )
)

server <- function(input, output, session) {
    # QC
    relative.abundances <- subset(vOTU.relative.abundance, select=-ID)
    sel.sample <- reactive({input$sel.sample})
    output$abundance.plot <- renderPlot({
        if ("All" %in% sel.sample()) {
            melt(relative.abundances) %>%
            ggplot(aes(x = variable, y=value)) +
            geom_boxplot(aes(fill=variable)) +
            labs(title = "Relative abundances", x = "Sample", y = "Relative abundance", fill="Sample") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- relative.abundances %>% select(sel.sample())
            melt(df) %>%
            ggplot(aes(x = variable, y=value)) +
            geom_boxplot(aes(fill=variable)) +
            labs(title = "Relative abundances", x = "", y = "Relative abundance", fill="Sample") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    quality.counts <- table(vOTU.stats.combined$checkv_quality)
    sel.quality <- reactive({input$sel.quality})
    output$checkv.quality.plot <- renderPlot({
        if ("All" %in% sel.quality()) {
            data.frame(quality.counts) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "CheckV quality", x = "Quality", y = "Number of vOTUs", fill="Quality") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- data.frame(quality.counts) %>% filter(Var1 %in% sel.quality())
            data.frame(df) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "CheckV quality", x = "Quality", y = "Number of vOTUs", fill="Quality") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    # Assembly
    sel.sample.assembly1 <- reactive({input$sel.sample.assembly1})
    output$assembly1.plot <- renderPlot({
        if ("All" %in% sel.sample.assembly1()) {
            melt(combined.sample.stats) %>%
            ggplot(aes(x = variable, y=value, group=Sample)) +
            geom_line(aes(colour=Sample), lwd=2) +
            labs(title = "Assembly statistics", x = "Workflow interval", y = "Number of contigs") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- combined.sample.stats %>% filter(Sample %in% sel.sample.assembly1())
            melt(df) %>%
            ggplot(aes(x = variable, y=value, group=Sample)) +
            geom_line(aes(colour=Sample), lwd=2) +
            labs(title = "Assembly statistics", x = "Workflow interval", y = "# contigs") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    # Taxonomy
    family.counts <- table(vOTU.stats.combined$Family)
    sel.family <- reactive({input$sel.family})
    output$family.plot <- renderPlot({
        if ("All" %in% sel.family()) {
            data.frame(family.counts) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified families", x = "Family", y = "Number of vOTUs", fill="Family") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- data.frame(family.counts) %>% filter(Var1 %in% sel.family())
            data.frame(df) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified families", x = "Family", y = "Number of vOTUs", fill="Family") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    subfamily.counts <- table(vOTU.stats.combined$Subfamily)
    sel.subfamily <- reactive({input$sel.subfamily})
    output$subfamily.plot <- renderPlot({
        if ("All" %in% sel.subfamily()) {
            data.frame(subfamily.counts) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified subfamilies", x = "Subfamily", y = "Number of vOTUs", fill="Subfamily") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- data.frame(subfamily.counts) %>% filter(Var1 %in% sel.subfamily())
            data.frame(df) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified subfamilies", x = "Subfamily", y = "Number of vOTUs", fill="Subfamily") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    genus.counts <- table(vOTU.stats.combined$Genus)
    sel.genus <- reactive({input$sel.genus})
    output$genus.plot <- renderPlot({
        if ("All" %in% sel.genus()) {
            data.frame(genus.counts) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified genera", x = "Genus", y = "Number of vOTUs", fill="Genus") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- data.frame(genus.counts) %>% filter(Var1 %in% sel.genus())
            data.frame(df) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified genera", x = "Genus", y = "Number of vOTUs", fill="Genus") +
            theme_bw(base_size = 20) +
            theme(axis.text.x=element_text(angle=45, vjust=0.9, hjust=1)) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    # Functional analysis
    type.counts <- table(vOTU.stats.combined$Type)
    sel.type <- reactive({input$sel.type})
    output$type.plot <- renderPlot({
        if ("All" %in% sel.type()) {
            data.frame(type.counts) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified types", x = "Type", y = "Number of vOTUs", fill="Type") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- data.frame(type.counts) %>% filter(Var1 %in% sel.type())
            data.frame(df) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified types", x = "Type", y = "Number of vOTUs", fill="Type") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    provirus.counts <- table(vOTU.stats.combined$provirus)
    sel.provirus <- reactive({input$sel.provirus})
    output$provirus.plot <- renderPlot({
        if ("All" %in% sel.provirus()) {
            data.frame(provirus.counts) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified proviruses", x = "Provirus", y = "Number of vOTUs", fill="Provirus") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        } else {
            df <- data.frame(provirus.counts) %>% filter(Var1 %in% sel.provirus())
            data.frame(df) %>%
            mutate(Var1 = str_replace_all(Var1, "n.a.", "Unclassified")) %>%
            ggplot(aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified proviruses", x = "Provirus", y = "Number of vOTUs", fill="Provirus") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2))
        }
    })
    # inStrain
    compared.samples$scaffold <- as.integer(str_replace(compared.samples$scaffold, "vOTU_", ""))
    compared.samples <- compared.samples %>% arrange(scaffold)
    compared.samples$name1 <- str_replace(compared.samples$name1, ".map.sorted.bam", "")
    compared.samples$name2 <- str_replace(compared.samples$name2, ".map.sorted.bam", "")
    compared.samples$vs <- paste(compared.samples$name1, compared.samples$name2, sep=" vs ")
    compared.samples.counts <- compared.samples %>% count(vs)
    compared.samples.means <- compared.samples %>%
    group_by(vs) %>%
    summarise_at(vars(coverage_overlap,
                        compared_bases_count,
                        percent_genome_compared,
                        length, consensus_SNPs,
                        population_SNPs),
                list(mean=mean, sd=sd))
    output$instrain.hist <- renderPlot({
        compared.samples %>%
        ggplot(aes(x=scaffold, fill=vs)) +
        geom_histogram() +
        geom_rug() +
        theme_bw(base_size = 20) +
        theme(plot.margin = margin(2,2,2,2, "cm"),
        axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2)) +
        labs(x="vOTU", y="Count") +
        guides(fill=guide_legend(title="Compared samples"))
    })
    output$mean.cov.overlap <- renderPlot({
        compared.samples.means %>%
        ggplot(aes(x=vs, y=coverage_overlap_mean)) +
            geom_col() +
            geom_errorbar(aes(ymin=coverage_overlap_mean - coverage_overlap_sd / sqrt(compared.samples.counts$n),
                            ymax=coverage_overlap_mean + coverage_overlap_sd / sqrt(compared.samples.counts$n)),
                        width=0.2) +
            geom_jitter(data=compared.samples, aes(x=vs, y=coverage_overlap),
                        width=0.1, alpha=0.3, col="blue") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2)) +
            labs(x="Compared samples", y="Mean coverage overlap") +
            coord_flip()
    })
    output$mean.bases.count <- renderPlot({
        compared.samples.means %>%
        ggplot(aes(x=vs, y=compared_bases_count_mean)) +
            geom_col() +
            geom_errorbar(aes(ymin=compared_bases_count_mean - compared_bases_count_sd / sqrt(compared.samples.counts$n),
                            ymax=compared_bases_count_mean + compared_bases_count_sd / sqrt(compared.samples.counts$n)),
                        width=0.2) +
            geom_jitter(data=compared.samples, aes(x=vs, y=compared_bases_count),
                        width=0.1, alpha=0.3, col="blue") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2)) +
            labs(x="Compared samples", y="Mean bases count") +
            coord_flip()
    })
    output$percent.genome.compared <- renderPlot({
        compared.samples.means %>%
        ggplot(aes(x=vs, y=percent_genome_compared_mean)) +
            geom_col() +
            geom_errorbar(aes(ymin=percent_genome_compared_mean - percent_genome_compared_sd / sqrt(compared.samples.counts$n),
                            ymax=percent_genome_compared_mean + percent_genome_compared_sd / sqrt(compared.samples.counts$n)),
                        width=0.2) +
            geom_jitter(data=compared.samples, aes(x=vs, y=percent_genome_compared),
                        width=0.1, alpha=0.3, col="blue") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2)) +
            labs(x="Compared samples", y="Mean % genome compared") +
            coord_flip()
    })
    output$mean.length <- renderPlot({
        compared.samples.means %>%
        ggplot(aes(x=vs, y=length_mean)) +
            geom_col() +
            geom_errorbar(aes(ymin=length_mean - length_sd / sqrt(compared.samples.counts$n),
                            ymax=length_mean + length_sd / sqrt(compared.samples.counts$n)),
                        width=0.2) +
            geom_jitter(data=compared.samples, aes(x=vs, y=length),
                        width=0.1, alpha=0.3, col="blue") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2)) +
            labs(x="Compared samples", y="Mean length") +
            coord_flip()
    })
    output$mean.consensus.SNPs <- renderPlot({
        compared.samples.means %>%
        ggplot(aes(x=vs, y=consensus_SNPs_mean)) +
            geom_col() +
            geom_errorbar(aes(ymin=consensus_SNPs_mean - consensus_SNPs_sd / sqrt(compared.samples.counts$n),
                            ymax=consensus_SNPs_mean + consensus_SNPs_sd / sqrt(compared.samples.counts$n)),
                        width=0.2) +
            geom_jitter(data=compared.samples, aes(x=vs, y=consensus_SNPs),
                        width=0.1, alpha=0.3, col="blue") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2)) +
            labs(x="Compared samples", y="Mean consensus SNPs") +
            coord_flip()
    })
    output$mean.pop.SNPs <- renderPlot({
        compared.samples.means %>%
        ggplot(aes(x=vs, y=population_SNPs_mean)) +
            geom_col() +
            geom_errorbar(aes(ymin=population_SNPs_mean - population_SNPs_sd / sqrt(compared.samples.counts$n),
                            ymax=population_SNPs_mean + population_SNPs_sd / sqrt(compared.samples.counts$n)),
                        width=0.2) +
            geom_jitter(data=compared.samples, aes(x=vs, y=population_SNPs),
                        width=0.1, alpha=0.3, col="blue") +
            theme_bw(base_size = 20) +
            theme(plot.margin = margin(2,2,2,2, "cm"), legend.position = "bottom",
            axis.title.x = element_text(vjust=-2), axis.title.y = element_text(vjust=2)) +
            labs(x="Compared samples", y="Mean population SNPs") +
            coord_flip()
    })
    # Stop app
    session$onSessionEnded(function() {
        stopApp()
    })
}

app <- shinyApp(ui = UI, server = server)
stopApp(app)
runApp(app, port=as.integer(port))
