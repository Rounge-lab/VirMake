options(warn=-1)

require("shiny")
require("shinydashboard")
require("shinythemes")
require("tidyverse")
require("ggplot2")
require("reshape2")


args = commandArgs(trailingOnly=TRUE)
stats.path <- args[1]
combined.sample.stats <- read.table(paste(stats.path, "Combined_Sample_stats.tsv", sep = "/"), header = TRUE, sep = "\t")
sample.stats.vibrant <- read.table(paste(stats.path, "Sample_stats_vibrant.tsv", sep = "/"), header = TRUE, sep = "\t")
sample.stats.virsorter2 <- read.table(paste(stats.path, "Sample_stats_virsorter2.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.AMGs <- read.table(paste(stats.path, "vOTU_AMGs.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.mapped.to.reads <- read.table(paste(stats.path, "vOTU_mapped_to_reads.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.relative.abundance <- read.table(paste(stats.path, "vOTU_Relative_Abundance.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.stats.combined <- read.table(paste(stats.path, "vOTU_stats_combined.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.stats.vibrant <- read.table(paste(stats.path, "vOTU_Stats_vibrant.tsv", sep = "/"), header = TRUE, sep = "\t")
vOTU.stats.virsorter2 <- read.table(paste(stats.path, "vOTU_Stats_virsorter2.tsv", sep = "/"), header = TRUE, sep = "\t")

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
                    selected = "All"
                ),
                checkboxGroupInput(
                    inputId = "sel.quality",
                    label = "Filter by quality:",
                    choices = sort(c("All", unique(vOTU.stats.combined$checkv_quality))),
                    selected = "All"
                ),
            ),
            mainPanel(
                plotOutput("abundance.plot"),
                plotOutput("checkv.quality.plot")
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
                plotOutput("assembly1.plot"),
                # plotOutput("assembly2.plot")
            )
        ),
        tabPanel("Taxonomy",
            titlePanel("Taxonomy results"),
            sidebarPanel(
                checkboxGroupInput(
                    inputId = "sel.family",
                    label = "Filter by family:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Family))),
                    selected = "All"
                ),
                checkboxGroupInput(
                    inputId = "sel.subfamily",
                    label = "Filter by subfamily:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Subfamily))),
                    selected = "All"
                ),
                checkboxGroupInput(
                    inputId = "sel.genus",
                    label = "Filter by genus:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Genus))),
                    selected = "All"
                ),
            ),
            mainPanel(
                plotOutput("family.plot"),
                plotOutput("subfamily.plot"),
                plotOutput("genus.plot"),
            )
        ),
        tabPanel("Functional analysis",
            titlePanel("Functional analysis results"),
            sidebarPanel(
                checkboxGroupInput(
                    inputId = "sel.type",
                    label = "Filter by type:",
                    choices = sort(c("All", unique(vOTU.stats.combined$Type))),
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
                plotOutput("type.plot"),
                plotOutput("provirus.plot"),
            )
        )
    )
)

server <- function(input, output) {
    # QC
    relative.abundances <- subset(vOTU.relative.abundance, select=-ID)
    sel.sample <- reactive({input$sel.sample})
    output$abundance.plot <- renderPlot({
        if ("All" %in% sel.sample()) {
            ggplot(melt(relative.abundances), aes(x = variable, y=value)) +
            geom_boxplot(aes(fill=variable)) +
            labs(title = "Relative abundances", x = "Sample", y = "Relative abundance", fill="Sample")
        } else {
            print(relative.abundances)
            df <- relative.abundances %>% select(sel.sample())
            ggplot(melt(df), aes(x = variable, y=value)) +
            geom_boxplot(aes(fill=variable)) +
            labs(title = "Relative abundances", x = "Sample", y = "Relative abundance", fill="Sample")
        }
    })
    quality.counts <- table(vOTU.stats.combined$checkv_quality)
    sel.quality <- reactive({input$sel.quality})
    output$checkv.quality.plot <- renderPlot({
        if ("All" %in% sel.quality()) {
            ggplot(data.frame(quality.counts), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "CheckV quality", x = "Quality", y = "Number of vOTUs", fill="Quality")
        } else {
            df <- data.frame(quality.counts) %>% filter(Var1 %in% sel.quality())
            ggplot(data.frame(df), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "CheckV quality", x = "Quality", y = "Number of vOTUs", fill="Quality")
        }
    })

    # Assembly
    sel.sample.assembly1 <- reactive({input$sel.sample.assembly1})
    output$assembly1.plot <- renderPlot({
        if ("All" %in% sel.sample.assembly1()) {
            ggplot(melt(combined.sample.stats), aes(x = variable, y=value, group=Sample)) +
            geom_line(aes(colour=Sample)) +
            labs(title = "Assembly statistics", x = "Workflow interval", y = "# contigs")
        } else {
            df <- combined.sample.stats %>% filter(Sample %in% sel.sample.assembly1())
            ggplot(melt(df), aes(x = variable, y=value, group=Sample)) +
            geom_line(aes(colour=Sample)) +
            labs(title = "Assembly statistics", x = "Workflow interval", y = "# contigs")
        }
    })

    # Taxonomy
    family.counts <- table(vOTU.stats.combined$Family)
    sel.family <- reactive({input$sel.family})
    output$family.plot <- renderPlot({
        if ("All" %in% sel.family()) {
            ggplot(data.frame(family.counts), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified families", x = "Family", y = "Number of vOTUs", fill="Family")
        } else {
            df <- data.frame(family.counts) %>% filter(Var1 %in% sel.family())
            ggplot(data.frame(df), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified families", x = "Family", y = "Number of vOTUs", fill="Family")
        }
    })
    subfamily.counts <- table(vOTU.stats.combined$Subfamily)
    sel.subfamily <- reactive({input$sel.subfamily})
    output$subfamily.plot <- renderPlot({
        if ("All" %in% sel.subfamily()) {
            ggplot(data.frame(subfamily.counts), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified subfamilies", x = "Subfamily", y = "Number of vOTUs", fill="Subfamily")
        } else {
            df <- data.frame(subfamily.counts) %>% filter(Var1 %in% sel.subfamily())
            ggplot(data.frame(df), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified subfamilies", x = "Subfamily", y = "Number of vOTUs", fill="Subfamily")
        }
    })
    genus.counts <- table(vOTU.stats.combined$Genus)
    sel.genus <- reactive({input$sel.genus})
    output$genus.plot <- renderPlot({
        if ("All" %in% sel.genus()) {
            ggplot(data.frame(genus.counts), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified genera", x = "Genus", y = "Number of vOTUs", fill="Genus")
        } else {
            df <- data.frame(genus.counts) %>% filter(Var1 %in% sel.genus())
            ggplot(data.frame(df), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified genera", x = "Genus", y = "Number of vOTUs", fill="Genus")
        }
    })

    # Functional analysis
    type.counts <- table(vOTU.stats.combined$Type)
    sel.type <- reactive({input$sel.type})
    output$type.plot <- renderPlot({
        if ("All" %in% sel.type()) {
            ggplot(data.frame(type.counts), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified types", x = "type", y = "Number of vOTUs", fill="Type")
        } else {
            df <- data.frame(type.counts) %>% filter(Var1 %in% sel.type())
            ggplot(data.frame(df), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified types", x = "type", y = "Number of vOTUs", fill="Type")
        }
    })
    provirus.counts <- table(vOTU.stats.combined$provirus)
    sel.provirus <- reactive({input$sel.provirus})
    output$provirus.plot <- renderPlot({
        if ("All" %in% sel.provirus()) {
            ggplot(data.frame(provirus.counts), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified proviruses", x = "Provirus", y = "Number of vOTUs", fill="Provirus")
        } else {
            df <- data.frame(provirus.counts) %>% filter(Var1 %in% sel.provirus())
            ggplot(data.frame(df), aes(x=Var1, y=Freq)) +
            geom_col(aes(fill=Var1)) +
            labs(title = "Identified proviruses", x = "Provirus", y = "Number of vOTUs", fill="Provirus")
        }
    })
}

app <- shinyApp(ui = UI, server = server)
runApp(app, port=8080)
