#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

using<-function(...) {
    libs<-unlist(list(...))
    req<-unlist(lapply(libs,require,character.only=TRUE))
    need<-libs[req==FALSE]
    if(length(need)>0){
        install.packages(need)
        lapply(need,require,character.only=TRUE)
    }
}
pckgs <- list("tidyverse", "pROC", "shiny", "patchwork", "ggbeeswarm", "cutpointr", "truncnorm")
#using(pckgs)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(cutpointr)
library(truncnorm)


# Define UI for application
ui <- fluidPage(responsive = FALSE,

    # Application title
    titlePanel("Biomarker metrics"),

    # Sidebar with a sliders
    sidebarLayout(
        sidebarPanel(
            sliderInput("npts",
                        "Number of subjects per group:",
                        min = 20,
                        max = 100,
                        value = 75),
            sliderInput("mean_diff",
                        "Difference in biomarker levels between healthy people and those with cancer:",
                        min = 0.5,
                        max = 3,
                        value = 1.75),
            sliderInput("popsd",
                        "Biological variation:",
                        min = 0.5,
                        max = 2,
                        value = 1),
            sliderInput("mnoise",
                        "Measurement noise [%]:",
                        min = 0,
                        max = 50.0,
                        value = 0),
            selectInput("mMetric", "Constrain metric:",
                        c("Max sensitivity given specificity" = "sens",
                          "Max specificity given sensitivity" = "spec"),
                        selected = "spec"),
            sliderInput("constr",
                        "Minimal value at [%]:",
                        min = 0,
                        max = 100,
                        value = 95),
            checkboxInput("checkCutoff", "Show cut-off value for detecting cancer", FALSE),
            checkboxInput("useBM", "Use biomarker", FALSE)
            #checkboxInput("checkShade", "Shade diagnosed as cancer", FALSE)

        ),

        # Show a plot of the generated distribution
        mainPanel(
            markdown("
    ### Detecting cancer using a biomarker
    Scenario: You have discovered a protein that is elevated in plasma of cancer patients and want to use it as a diagnostic biomarker.
    Explore different factors:

    * Biological differences of interest (cancer vs healthy)
    * Biological variability (how biomarker levels vary in a population)
    * Measurement/technical noise (how well your biomarker measurement works)
    * Discrimination threshold (decision on biomarker level that indicates cancer)

    ROC curve shows you the overall discriminatory performance of the biomarker.
    The table below shows other metrics (AUC, sensitivity, specificity, accuracy).

    Select '**Use biomarker**' to see how well it would perform in a clinical setting.
    "),

            plotOutput("distPlot"),
            tableOutput("myTable"),
            textOutput("myTableCMC"),
            textOutput("myTableCMH")

        )
    )
)

server <- function(input, output) {
    create_data <- reactive({
        set.seed(100)
        npts <- input$npts
        mean_diff <- input$mean_diff
        popsd <- input$popsd
        mnoise <- input$mnoise
        cdata <- data.frame(
            group = factor(rep(c("Yes","No"), each = npts), levels = c("No","Yes")),
            biomarker = c(rtruncnorm(npts, a = 0, mean = 1.5 + mean_diff, sd = popsd),
                          rtruncnorm(npts, a = 0, mean = 1.5, sd = popsd)),
            measnoise = c(rnorm(2*npts, mean = 1, sd = mnoise/100))
        ) %>% mutate(biomarker = biomarker * measnoise)

    })

    calculate_metrics <- reactive({
        cdata <- create_data()
        mconstrain <- input$constr/100
        mMetric <- input$mMetric
        x <- cdata$biomarker
        class <- cdata$group
        cp <- cutpointr(x = x,
                        class = class,
                        pos_class = "Yes", neg_class = "No", direction = ">=",
                        method = maximize_metric,
                        metric = ifelse(mMetric == "sens", sens_constrain, spec_constrain),
                        min_constrain = mconstrain)
        return(cp)
    })

    output$distPlot <- renderPlot({
        cdata <- create_data()
        cmet <- calculate_metrics()

        thr <- summary(cmet)$cutpointr[[1]]$optimal_cutpoint
        thrSpec <- summary(cmet)$cutpointr[[1]]$specificity
        thrSens <- summary(cmet)$cutpointr[[1]]$sensitivity
        p0 <- cdata %>% ggplot(aes(x = 1, y = biomarker, color = group)) +
            theme_bw()+
            ggbeeswarm::geom_beeswarm(cex =2)+

            ylab("Biomarker value") +
            scale_color_brewer(palette = "Set1", direction = -1)+
            theme(legend.position = "none") + scale_x_discrete() +xlab("")

        p1 <- cdata %>% ggplot(aes(x = group, y = biomarker, color = group)) +
            theme_bw()+
            ggbeeswarm::geom_beeswarm(cex =2)+
            xlab("Disease") +
            ylab("Biomarker value") +
            scale_color_brewer(palette = "Set1", direction = -1)+
            theme(legend.position = "none", aspect.ratio = 1)


        croc <- pROC::roc(group~biomarker, data = cdata)
        p2<- pROC::ggroc(croc, size = 2) + ggtitle("ROC curve") +
            geom_abline(slope = 1, intercept = 1, linetype = "dashed") +
            scale_y_continuous(expand = c(0.001,0.001)) + scale_x_reverse(expand = c(0.001,0.001))+
            theme_bw() + coord_equal() +
            xlab("Specificity") + ylab("Sensitivity")


        if (input$checkCutoff){
            p1 <- p1 + geom_hline(yintercept = thr,
                                  linetype = "dashed")
            p2 <- p2 + annotate(geom = "point",
                                x = thrSpec,
                                y = thrSens,
                                size = 5,
                                shape = 21,
                                stroke = 3,
                                color = "orange")
        }

        p3 <- p1 + p2
        p3

    })

    output$myTable = renderTable( align = "c",{
        cp <- summary(calculate_metrics())
        cp_tab<- cp$cutpointr[[1]] %>% mutate(
            Specificity = 100*specificity,
            Sensitivity = 100*sensitivity,
            Accuracy = 100*acc) %>%
            select(AUC, `Cut-off` = optimal_cutpoint, `Specificity %`=Specificity, `Sensitivity %`=Sensitivity, `Accuracy %`=Accuracy)
        print(cp_tab)
    })

    output$myTableCMC = renderText({
        cp <- summary(calculate_metrics())$confusion_matrix[[1]][,-1]
        if (input$useBM){
            sprintf("Out of %i subjects with cancer, you detected %i and %i remain undetected.", input$npts, cp$tp, cp$fn)
        }else sprintf("")
    })

    output$myTableCMH = renderText({
        cp <- summary(calculate_metrics())$confusion_matrix[[1]][,-1]
        if (input$useBM){
            sprintf("You told %i (of %i) healthy people they have cancer.", cp$fp, input$npts)
        } else sprintf("")
    })
}

# Run the application
shinyApp(ui = ui, server = server)
