require("shiny")
require("shinydashboard")
require("shinythemes")

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
            titlePanel("Assembly results")
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
        tabPanel("Functional analysis"),
    )

)

server <- function(input, output) {

    # QC
    relative.abundances <- subset(vOTU.relative.abundance, select=-ID)
    sel.sample <- reactive({input$sel.sample})
    output$abundance.plot <- renderPlot({
        if ("All" %in% sel.sample()) {
            boxplot(relative.abundances, main = "Relative abundances", xlab = "Sample", ylab = "Relative abundance")
        } else {
            boxplot(relative.abundances[,sel.sample()], main = "Relative abundances", xlab = "Sample", ylab = "Relative abundance")
        }
    })
    quality.counts <- table(vOTU.stats.combined$checkv_quality)
    sel.quality <- reactive({input$sel.quality})
    output$checkv.quality.plot <- renderPlot({
        if ("All" %in% sel.quality()) {
            barplot(quality.counts, main = "CheckV quality", xlab = "Quality", ylab = "Number of vOTUs")
        } else {
            barplot(quality.counts[sel.quality()], main = "CheckV quality", xlab = "Quality", ylab = "Number of vOTUs")
        }
    })

    # Taxonomy
    family.counts <- table(vOTU.stats.combined$Family)
    sel.family <- reactive({input$sel.family})
    output$family.plot <- renderPlot({
        if ("All" %in% sel.family()) {
            barplot(family.counts, main = "Identified families", xlab = "Family", ylab = "Number of vOTUs")
        } else {
            barplot(family.counts[sel.family()], main = "Identified families", xlab = "Family", ylab = "Number of vOTUs")
        }
    })
    subfamily.counts <- table(vOTU.stats.combined$Subfamily)
    sel.subfamily <- reactive({input$sel.subfamily})
    output$subfamily.plot <- renderPlot({
        if ("All" %in% sel.subfamily()) {
            barplot(subfamily.counts, main = "Identified subfamilies", xlab = "Subfamily", ylab = "Number of vOTUs")
        } else {
            barplot(subfamily.counts[sel.subfamily()], main = "Identified subfamilies", xlab = "Subfamily", ylab = "Number of vOTUs")
        }
    })
    genus.counts <- table(vOTU.stats.combined$Genus)
    sel.genus <- reactive({input$sel.genus})
    output$genus.plot <- renderPlot({
        if ("All" %in% sel.genus()) {
            barplot(genus.counts, main = "Identified genera", xlab = "Genus", ylab = "Number of vOTUs")
        } else {
            barplot(genus.counts[sel.genus()], main = "CheckV quality", xlab = "Genus", ylab = "Number of vOTUs")
        }
    })


}

shinyApp(ui = UI, server = server)
