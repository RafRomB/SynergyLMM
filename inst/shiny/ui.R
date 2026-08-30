# Define UI
ui <- fluidPage(
  shinyjs::useShinyjs(), # Include shinyjs
  theme = bslib::bs_theme(bootswatch = "united", bg = "white", fg = "#171d1f",
                          primary = "#29454d", secondary = "#2c778a"),
  
  tags$style(HTML("
    /* General styles */
    body {
      font-family: 'Ubuntu', sans-serif;
    }

    /* Style for the title panel */
    .title-panel {
      background-color: #a7c9c7; /* Background color */
      color: #ffffff; /* Text color */
      padding: 10px; /* Add some padding */
      border-radius: 10px; /* Optional: Rounded corners */
      font-size: 50px;
      text-align: left;
    }

    /* Style for the top tab bar */
    .nav-tabs {
      background-color: #a7c9c7; /* Background color */
      border-radius: 10px; /* Optional: Rounded corners */
    }

    /* Style for the active tab */
    .nav-tabs .nav-item .nav-link.active {
      background-color: #2c778a; /* Active tab background */
      color: white; /* Active tab text color */
      font-weight: bold; /* Bold font for active tab */
    }

    /* Style for non-active tabs */
    .nav-tabs .nav-item .nav-link {
      color: #f5f5f5; /* Non-active tab text color */
      background-color: #3b6475; /* Non-active tab background color */
      font-weight: bold;
    }

    /* Hover effect for non-active tabs */
    .nav-tabs .nav-item .nav-link:hover {
      background-color: #1f2a30; /* Hover background color */
      color: #ffffff; /* Hover text color */
    }

    /* Floating sidebar styles for larger screens */
    .sidebar-column {
      background-color: #e8f4f3;
      padding: 15px;
      border-radius: 8px;
    }
    
    /* --- stick the sidebar on large screens --- */
    .floating-sidebar {
      position: -webkit-sticky;
      position: sticky;
      top: 80px;
      align-self: flex-start;
      max-height: calc(100vh - 100px);
      overflow-y: auto;
    }

    .main-panel-adjust {
      margin-top: 20px;
    }

    /* reset on small screens */
    @media (max-width: 1200px) {
      .floating-sidebar {
        position: static;
        max-height: none;
        overflow: visible;
      }

      .main-panel-adjust {
        margin-left: 0;
      }

      .nav-tabs {
        font-size: 14px;
      }

      .title-panel {
        font-size: 22px;
      }
    }

    @media (max-width: 768px) {
      .floating-sidebar {
        position: static;
        width: 100%;
        margin-bottom: 10px;
      }

      .main-panel-adjust {
        margin-left: 0;
      }

      .title-panel {
        font-size: 20px;
        text-align: left;
      }

      .nav-tabs .nav-item .nav-link {
        font-size: 12px;
      }
    }

    @media (max-width: 576px) {
      .sidebar-column {
        padding: 10px;
        font-size: 14px;
      }

      .title-panel {
        font-size: 18px;
        padding: 8px;
      }

      .main-panel-adjust {
        margin-top: 10px;
      }
    }
")),
  
  tags$div(class = "title-panel",
           tags$img(
             src = "SynergyLMM_icon.png",   # Put your icon file in the www folder
             height = "75px",    # Adjust size
             style = "vertical-align: middle; margin-right: 15px;"
           ),
           tags$span(strong(paste0("SynergyLMM (v", utils::packageVersion("SynergyLMM"), ")")))),
  
  #Tab layout
  tabsetPanel(
    
    # 0. Intro Tab for method description ----
    tabPanel("Tutorial",
             sidebarLayout(
               sidebarPanel(
                 tags$h4("Tutorial Sections"),
                 tags$ol(
                   tags$li(tags$a(href = "#overview", "SynergyLMM Overview")),
                   tags$li(tags$a(href = "#upload", "Data Upload and Model Estimation")),
                   tags$li(tags$a(href = "#synergy", "Synergy Analysis")),
                   tags$li(tags$a(href = "#diagnostics", "Model Diagnostics")),
                   tags$ul(
                     tags$li(tags$a(href = "#influential", "Influential Diagnostics"))
                     ),
                   tags$li(tags$a(href = "#performance", "Model Performance")),
                   tags$li(tags$a(href = "#power_analysis", "Power Analysis")),
                   tags$ul(
                     tags$li(tags$a(href = "#posthoc", "Post Hoc Power Analysis")),
                     tags$li(tags$a(href = "#apriori", "A Priori Power Analysis")),
                     tags$ul(
                       tags$li(tags$a(href = "#samplesize", "Sample Size")),
                       tags$li(tags$a(href = "#timepwr", "Time and Frequency")),
                       tags$li(tags$a(href = "#variability", "Data Variability"))
                       )
                     ),
                   tags$li(tags$a(href = "#workflow", "Recommended workflow for using SynergyLMM"))
                   ), 
                 width = 3, 
                 class = "floating-sidebar",
               ), 
               mainPanel(
                 class = "main-panel-adjust",
                 card(id = "overview",
                      card_header(h2("SynergyLMM Workflow Overview")),
                 tags$figure(
                   tags$img(
                     width = "100%",
                     src = "LMM_Graphical_Abstract.png",
                     alt = "SynergyLMM workflow overview."
                 ),
                 tags$figcaption(tags$b("Overview of SynergyLMM workflow. a, "),
                                 "The first step is uploading the longitudinal tumor burden-related measurements for the different treatment groups.
                                 The input consists of a tabular data in long format with at least 4 columns containing information abouth the samples IDs, the time points for each measurement, the treatment group, and the tumor measurement.",
                                 tags$b("b, "),
                                 "The input data will then be processed to estimate the (non-)linear mixed effect model that fits the tumor growth dynamics for each subject, and estimate the growth rates for each treatment group.",
                                 tags$b("c, "), 
                                 "SynergyLMM offer various functions to evaluate the model diagnostics and model performance, as well as for the identification of potential outliers and influential individuals for the model fit and the treatment groups.",
                                 tags$b("d, "),
                                 "Once a proper model with satisfactory diagnostics has been obtained, the statistical assessment of combination effects is performed, with a time dependent estimation of the synergy score and combination index, along with their confidence intervals and statistical significance. 
                                 The method allows for testing synergy using three reference models: Bliss independence, highest single agent (HSA), and response additivity (RA).",
                                 tags$b("e, "),
                                 "SynergyLMM implements versatile functions to calculate the", 
                                 tags$i("post hoc"),
                                 "power of the experiment and the",
                                 tags$i("a priori"), 
                                 "power by modifying the values of key experimental variables, such as sample size,  follow-up time, or frequency of the measurements.")
                 
                 )
               ),
               card(id = "example_data",
                 card_header(h2("How to Use SynergyLMM")),
                 tags$p(tags$b("Welcome to the SynergyLMM tutorial!"), tags$br(),
                        "This quick-start tutorial shows how to use the web-app.", "An example dataset will be used,",
                        "which can be downloaded by clicking in the link below:"),
                 tags$a(href = "eg_grwth_data.xlsx", download = NA, "Download Example Dataset")
                 
               ),
               card(id = "upload",
                    card_header(h3("Data Upload and Model Fit")),
                    tags$p("The user can upload the data in the", tags$b("'Model Estimation'"), "tab.",
                           "The input data should be a data table in long format with at least four columns: one for sample IDs,",
                           "one for time points, one for treatment groups, and one for tumor volume measurements.",
                           "(See the example dataset for reference)."), 
                    tags$p("This is how the data should look like:"),
                    card(h3("Input Table"), DTOutput("ex_dt")),
                    tags$p("Once the data is uploaded, the user will need to",
                           "specify which columns contain each information. For this example, columns", tags$i("Mouse, Day, Group,"),
                           "and", tags$i("TumVol"), "contain the information about the sample IDs, time points, treatments, and tumor measurements,",
                           "respectively (Fig. 1a)."),
                    tags$p("The next step is to define the time points for the experiment. By default, the minimum and maximum values",
                           "in the 'Time' column will be used as the start and end time points, but these can be modified.",
                           "The name of the different treatments in the experiment also need to be specified. In this example,",
                           "'Control', 'DrugA', 'DrugB', and 'Combination' correspond to the control, single drugs, and the drug",
                           "combination groups (Fig. 1b)."),
                    tags$p("Next, clicking on 'Run Analysis' will fit the model and produce a table with the model estimates,",
                           "(Fig. 1c), including the", tags$i("tumor growth rates"), " of the different groups, and the standard deviations",
                           "for the random effects and for the residuals.", "If 'Show Plot' is selected, the plots with the tumor growth curves",
                           "are displayed along with the regression line for the fixed effect coefficients (i.e., the tumor growth rates for the control and treatment groups, Fig. 1d)."),
                    tags$p("The orange-shaded lines in Fig. 1d indicate two individuals, mouse 2 and mouse 10, whose measurements are (intentionally) notably different from the others,",
                           "and which will greatly affect the model and the results. The impact of these individuals and how to deal with them will be discussed later."),
                    tags$figure(
                      tags$img(
                        width = "100%",
                        src = "Fig1.png"
                      ),
                      tags$figcaption(tags$b("Figure 1. a,"), 
                                             "Upload dataset.",
                                             tags$b("b, "), 
                                             "Set time and name of treatments.",
                                             tags$b("c, "),
                                             "Table with model estimates.",
                                             tags$b("d, "),
                                             "Plots of tumor growth data. The arrow heads indicate subjects for which only four observations were made. 
                                             The orange-shaded lines indicate individuals whose measurements are notably different from the others.")),
                    tags$p(tags$b("Minimum Observations per Subject and Growth Model"), 
                           tags$br(), tags$br(),
                           "In Fig. 1b, there is an option called", tags$b("'Minimum Observations per Subject'"),", which allows for excluding",
                           "samples with fewer than the specified number of observations. For example, in Figure 1d, the arrow heads indicate",
                           "subjects for which only four observations were made. If the the minimum observations per subject is set to five,",
                           "these samples will be omitted from the analysis (Fig 2).", 
                           tags$br(), tags$br(),
                           "Fig. 1b also shows an option called", tags$b("'Growth Model'"),". This box allows to choose between 'Exponential' or 'Gompertz' growth models. The exponential
                           model is fitted using a linear mixed effect model, while the Gompertz model is fitted using a non-linear mixed effect model.
                           The exponential model is generally a good starting point for most in vivo studies, and it tends to converge more easily. However, it may be 
                           insufficient when tumor growth dynamics deviate significantly from exponential behavior. In such cases, the Gompertz model typically 
                           offers better flexibility and more accurate estimates of tumor growth and treatment effects."),
                    tags$figure(
                      tags$img(
                        width = "100%",
                        src = "Fig2.png"
                      ),
                      tags$figcaption(tags$b("Figure 2. a, "), 
                                      "Choosing subjects with at least five observations.",
                                      tags$b("b, "), 
                                      "Plots of tumor growth data from which the individuals with less than five observations have been removed.",
                                      tags$b("c, "),
                                      "Choosing the Gompertz growth model.",
                                      tags$b("d, "), 
                                      "Plots of tumor growth data using the Gompertz growth model. Note the difference in the regression lines in the plots shown in", tags$b("a"), ".")
                      )
                    ),
               
               card(id = "synergy",
                    card_header(h3("Synergy Analysis")),
                    tags$p("Once the model is fitted, the",
                           tags$b("Synergy Analysis"), "can be performed.",
                           "Fig. 3a shows the different options for synergy analysis.",
                           "The user can choose between Bliss, highest single agent (HSA), or response additivity (RA)",
                           "reference models. In the advance options, selecting 'Use Robust Estimators' will apply",
                           "sandwich-based cluster-robust variance estimators. This is recommended and can help to deal with",
                           "possible model misspecifications."),
                    tags$p("The analysis for the response additivity reference model is based on simulations, and the number of",
                           "simulations to run can also be specified by the user."),
                    tags$p("Fig. 3a also shows two additional options for the statistical assessment of synergy, 'Confidence level' and 'Correction method for adjusting p-values'.",
                           "The first one allows the user to choose the confidence level for the confidence intervals and p-values calculation (default to 0.95),",
                           "while the second one allows the user to choose between different methods for multiple comparison synergy p-values adjustment.",
                           "More information about the available methods can be obtained by clicking in the '?' help menu."),
                    tags$p("The output of the analysis includes a plot with the combination index (CI) and synergy score (SS) for each",
                           "time point (Fig. 3b), as well as a table with all the results. The plot in Fig. 3b shows the estimated synergy metric (CI or SS) value,
                           together with the confidence interval (vertical lines), and the color of the dots indicates the statistical significance."),
                    tags$p("The interpretation of the CI and SS values follows the typical definition of CI and SS:"),
                    tags$ul(
                      tags$li(tags$b("CI:"),"Provides information about the observed drug combination effect versus the expected", 
                              "additive effect according to the reference synergy model. A drug combination effect larger",
                              "than the expected (CI > 1) indicates synergism, a drug combination effect equal to the",
                              "expected (CI =1) indicates additivity, and a lower drug combination effect than the expected",
                              "(CI < 1) indicates antagonism.", "The value of the CI represents the proportion of tumor cell",
                              "survival in the drug combination group compared to the expected tumor cell survival according",
                              "to the reference model."),
                      tags$li(tags$b("SS:"),"The synergy score is defined as the excess response due to drug interaction compared to",
                              "the reference model. Following this definition, a SS>0, SS=0, and SS<0, represent synergistic, additive",
                              "and antagonistic effects, respectively.")
                    ),
                    tags$figure(
                      tags$img(
                        width = "100%",
                        src = "Fig3.png"
                      ),
                      tags$figcaption(tags$b("Figure 3. a, "),
                                      "Side panel for running the synergy analysis.",
                                      tags$b("b, "),
                                      "Plot of synergy results. Each dot represents
                                      the point estimate for the Combination Index (CI) and Synergy Score (SS), while the vertical
                                      gray lines represent the 95% confidence interval. The color of the dots reflects the minus log10 p-value for
                                      the synergy test. The horizontal dashed lines indicate the thresholds for assessing synergy, corresponding to 1
                                      for the CI and 0 for the SS.")
                    ),
                    card(h3("Synergy Results"), DTOutput("ex_syn")),
                    tags$p("As observed in Fig. 3b and the table with the results, the only time point for which there is a statistically significant",
                           "result is time 3. At this time point, the CI is lower than 1, and the SS is higher than 0, indicating",
                           "synergy in the drug combination group. However, at the other time points, the results show that the confidence",
                           "interval includes 1 for the CI and 0 for the SS, and the p-value is higher than 0.05, indicating that",
                           "the additive effect of the drug combination cannot be rejected.", 
                           tags$br(), tags$br(),
                           "In addition to the synergy results, a table output with the model estimates is also returned. In the case of using the exponential model, 
                           as a model is fitted for each time window being analyzed, the model estimates for each time point are reported:"),
                    card(h3("Model Estimates by Time"), DTOutput("ex_syn_est"))
                    
               ),
               card(id = "diagnostics",
                    card_header(h3("Model Diagnostics")),
                    tags$p("An important step in the analysis is to verify that the model is appropriate for",
                           "the input data. In the", tags$b("Model Diagnostics"), " tab, the user can check the main",
                           "assumptions: the normality of the random effects and the residuals."), 
                    tags$p("The model diagnostics output includes four plots to help",
                           "evaluate the adequacy of the model. The two plots at the top of Fig. 4a show the Q-Q plot of the random effects and normalized residuals.",
                           "Most points in these plots should fall along the diagonal line.",
                           "The two plots at the bottom of Fig. 4a show scatter plots of the standardized residuals versus fitted values and standardized residuals per
                           time and per treatment, which give information about variability of the residuals and possible outlier observations.",
                           "The residuals in these plots should appear randomly distributed around the horizontal line, indicating homoscedasticity.",
                           "Fig. 4b shows the output of Shapiro-Wilk normality tests for the random effects and the residuals.
                           A non-significant p-value obtained by the Shapiro-Wilk normality indicates that there is no evidence to reject normalilty."), 
                    tags$figure(
                      tags$img(
                        width = "100%",
                        src = "Fig4.png"
                      ),
                      tags$figcaption(tags$b("Figure 4. a, "), 
                                      "Q-Q plots of the",
                                      "random effects (top left), normalized residuals (top right), 
                                      scatter plots of normalized residuals versus fitted values (bottom left), and normalized residuals per time and treatment group (bottom right).", 
                                      tags$b("b, "), 
                                      "Results of Shapiro-Wilk's normality test for the random effects and normalized residuals.")),
                    tags$p("The results from the plots and the normality tests seem to indicate some departure from the normality, especially evident for the residuals.",
                           "This could be related to those individuals which has a large variation in the measurements."),
                    tags$p("Observations with absolute normalilzed residuals greater", 
                           "than the 1−0.05/2 quantile of the standard normal distribution are identified as potential outlier observations. A table",
                           "with these observations is also provided by SynergyLMM. It can be observed that many observations",
                           "from mouse 2 and mouse 10 have large residual values."),
                    card(h4("Potential Outlier Observations"), DTOutput("ex_outlier")),
                    tags$p("As mentioned before, mouse 2 and mouse 10 are intentionally introduced outliers in the model, and there are many observations from these subjectes that are
                           indeed identified as potential outliers.",
                           "More details about outlier identification are provided in the",
                           tags$b("Influential Diagnostics"), "section. But first, there are several approaches that the user can follow to try to improve the model diagnostics."),
                    tags$p("If the diagnostic plots and tests show evident violations of the model assumptions, there are several solutions that may help improving the model:",
                           tags$ul(
                             tags$b("- Transform the time variable: "),
                             "Both exponential and Gompertz models are based on ordinary differential equations (ODEs),
                                   which can be formulated using arbitrary time scales. Applying transformations, such as the square root or logarithm (with +1 offset),
                                   can often improve model fit, enhance linearity, stabilize variance, and ensure better numerical performance. Importantly, when comparing
                                   the exact values of a drug combination synergy index, the models should use the same timescale"
                           ), tags$ul(
                             tags$b("- Specify a residual variance structure: "),
                             "SynergyLMM uses the ‘nlme’ R package, which supports flexible variance and correlation structures.
                                   These allow for heterogeneous variances or within-group correlations, such as those based on subject, time point, treatment, or combinations of these.
                                   Applying these structures can substantially improve model diagnostics and the robustness of the estimates. This can be done in the ",
                             
                           tags$i("Advanced Options"), "of the",
                           tags$b("Model Estimation"), "tab."),
                           tags$ul(tags$b("- Prefer the Gompertz model when the exponential model fails: "), 
                           "While the exponential model is generally a good starting point for most 
                                   in vivo studies, and it tends to converge more easily, it may be insufficient when tumor growth dynamics deviate significantly from exponential 
                                   behavior. In such cases, the Gompertz model typically offers better flexibility and more accurate estimates of tumor growth and treatment effects."),
                           tags$ul(tags$b("- Carefully address potential outlier: "), 
                           "Individuals or measurements highlited as potential outliers may warrant 
                                   further investigation to reveal the reasons behind unusual growth behaviours, and potentially exclude these before re-analysis, 
                                   after careful reporting and justification.")
                           ),
                    tags$p("In this case, the model will be fitted especifying a unequal variance for the errors. Concretely, a different variance for each treatment is defined (Fig. 5a).",
                           "It can be seen that this has improved the model diagnostics, and now both the random effects and normalized residuals are approximately normally distributed.",
                           "There are still some potential outlier observations, which is normal, since there are two intentional outlier subjects."),
                    tags$figure(
                      tags$img(
                        width = "100%",
                        src = "Fig5.png"
                      ),
                      tags$figcaption(tags$b("Figure 5. a, "),
                                      "Using the Advanced Options to specify a treatment-specific residual variance structure",
                                      tags$b("b, "),
                                      "Q-Q plots of the",
                                      "random effects (top left), normalized residuals (top right), 
                                      scatter plots of normalized residuals versus fitted values (bottom left), and normalized residuals per time and treatment group (bottom right) of the corrected model.", 
                                      tags$b("c, "), 
                                      "Results of Shapiro-Wilk's normality test for the random effects and normalized residuals.")),
                    card(id = "influential",
                         card_header(h3("Influential Diagnostics")),
                         tags$p(
                           "The",
                           tags$b("Influential Diagnostics"),
                           "tab is useful to identify",
                           "subjects that have a great impact in the model and that may warrant more careful analysis.",
                           "To this end, SynergyLMM uses subject-specific Cook's distances."
                         ), tags$p(
                           "Cook's distances can be calculated based on the change in fitted values, or based on the change in the fixed effects."
                         ), tags$p(
                           "Cook's distances based on the change in fitted values assess the influence of each subject on the model fit.
                                On the other hand, Cook's distance based on the change in fixed effect values assess the influence of each subject
                                in the estimated parameters of the treatment groups."
                         ), tags$p(
                           "Fig. 6 shows results of the influential diagnostics of our example model. It can be seen how Mouse 10 is identified as having a great influence in the fitted values."
                         ), tags$figure(
                           tags$img(
                             width = "100%",
                             src = "Fig6.png"
                           ),
                           tags$figcaption(tags$b("Figure 6. a, "),
                                           "Advanced Options for Model Diagnostics. The users can choose between different normality test, the threshold for outlier identification and
                                           for Cook's distances, and which type of Cook's distance to calculate.",
                                           tags$b("b, "),
                                           "Plot of the Cook's distance based on the change in fitted values versus subjects.", 
                                           "The highlighted points indicate those individuals whose Cook's distance is higher than three time the mean of the values."))
                         ),
                    ),
               card(id = "performance",
                    card_header(h3("Model Performance")),
                    tags$p("Another important tool for the diagnostics of the model is the plot of the observed vs. predicted values, available",
                           "in the", tags$b("Model Performance"),"tab.", 
                           "The output from this tab includes a table with model performance metrics, such as R-square",
                           "and various error measures (Fig. 7a), along with the plot of observed vs. predicted values (Fig. 7b)."),
                    tags$p("Each plot in Fig. 7b corresponds to a subject. The blue dots represent the actual measurements, the continuous line",
                           "shows the group (treatment) marginal predicted values, and the dashed line indicates each individual's",
                           "predicted values"),
                    tags$p("By examining how closely the lines fit the data points, the user can assess how well the model fits",
                           "the data.", "In this example, the model generally fits the data well, except for certain subjects, such as",
                           "mouse 2 and mouse 10."),
                    tags$figure(
                      tags$img(
                        width = "100%",
                        src = "Fig7.png"
                      ),
                      tags$figcaption(tags$b("Figure 7. a, "), "Table with some metrics about the model performance.", 
                                      tags$b("b, "), 
                                      "Plots of observed vs. predicted values for each subject. The dots indicate the actual measurements, 
                                      while the continuous and dashed lines show the marginal and subject-specific predicted values."))
                    ),
               
               card(id = "power_analysis", h2("Power Analysis"),
                    tags$h5(tags$b("Note: Power analyses are only available for models fitted using the exponential growth model.")), 
                    card(id = "posthoc",
                         card_header(h3("Post Hoc Power Analysis")),
                         tags$p("The", tags$b("Post Hoc Power Analysis"), "allows users to check the statistical power of the dataset and fitted model",
                                "for synergy analysis using the Bliss and HSA reference models."),
                         tags$p("The post hoc power calculation is based on simulations, and depending on the dataset, it can take",
                                "some time to run. The power is returned as the proportion of simulations resulting in a significant synergy",
                                "hypothesis test."),
                         tags$p("Fig. 8a shows the side panel for running the post hoc power analysis. The user can define the number of",
                                "simulations to run, the p-value threshold for considering significance, and the time point for which to",
                                "perform the analysis."),
                         
                         tags$figure(
                           tags$img(
                             width = "100%",
                             src = "Fig8.png"
                           ),
                           tags$figcaption(tags$b("Figure 8. a, "),"Side panel for running the post hoc power analysis.", 
                                           tags$b("b, "),
                                           "Results of the post hoc power analysis.")
                         )
                    ),
                    card(id = "apriori",
                         card_header(h3("A Priori Power Analysis")),
                         tags$p("The user can also perform three different types of a priori power",
                                "analyses by modifying various experimental variables. To do this, 'SynergyLMM' creates exemplary datasets using the estimates from the fitted model",
                                "in the", tags$b("Model Estimation"), "tab.",
                                "Then, several parameters can be modified and evaluated by the user, such as: the number of subjects per group, the days of measurements,",
                                "the coefficients (tumor growth rates) for each treatment group, the standard deviation of the random effects (between-subject variance),",
                                "and the standard deviation of the residuals (within-subject variance)."),
                         tags$p("SynergyLMM then provides the a priori power calculation for each exemplary dataset defined by a specific set of experimental parameters."),
                         card(tags$h4(id = "samplesize", "Sample Size"),
                              tags$p("The first option that can be tested is how the power would vary as the sample size per group changes.",
                                     "To do this, an exemplary dataset is used to fit models with the same values for all parameters (coefficients for the different",
                                     "treatment groups, standard deviations of the random effects and residuals, number of time points, etc.) as in the fitted model,",
                                     "but varying the sample size in each group. The user can select the range of sample sizes to evaluate and the reference model to use (Fig. 9a)."),
                              tags$p("When clicking 'Run Sample Size Power Analysis', a plot with the exemplary data and another showing the power variation with the",
                                     "sample size are displayed (Fig. 9b). It can be confirmed that the values shown in the exemplary data plot are the same as the",
                                     "model estimates shown in Fig. 1c."),
                              tags$figure(
                                tags$img(
                                  width = "100%",
                                  src = "Fig9.png"
                                ),
                                tags$figcaption(tags$b("Figure 9. a, "), 
                                                "Side panel for performing power analysis by varying the sample size per group.", 
                                                tags$b("b, "), 
                                                "Plot displaying how the statistical power varies 
                                                  with sample size per group. The dashed line indicates the statistical power threshold of 0.8.")
                                
                              ),
                              tags$p("The results from the a priori power analysis show that to achieve a statistical power of 0.8, and assuming the other parameters remain constant,",
                                     "the sample size per group should be at least 13.")),
                         card(tags$h4(id = "timepwr", "Time Power Analysis"),
                              tags$p("Another parameter that can affect statistical power is the duration of follow-up or the frequency of measurements.",
                                     "Given the results of a pilot experiment, the user may wonder how the results would change if tumor growth were monitored for",
                                     "a longer time period, or how often tumor growth measurements should be taken. To address this, the user can evaluate statistical",
                                     "power by varying the maximum follow-up time or the frequency of measurements, while keeping the other parameters fixed,",
                                     "as determined from the fitted model."),
                              tags$h5("Maximum Time of Follow-up"),
                              tags$p("First, it can be examined how power changes if the follow-up time is increased or decreased. To do this, select",
                                     "the analysis based on the maximum time of follow-up ('max'), as shown in Fig. 10a. The sample size",
                                     "per group, the time interval between measurements, and the reference model to use can be also specified. For this example, the analysis evaluated",
                                     "how power varied from 3 to 42 days of follow-up, with 5 subjects per group, and measurements",
                                     "taken every 3 days."),
                              
                              tags$figure(
                                tags$img(
                                  width = "100%",
                                  src = "Fig10.png"
                                ),
                                tags$figcaption(tags$b("Figure 10. a, "), 
                                                "Side panel for conducting power analysis by varying the maximum time of follow-up. Note that the type of analysis is set to 'max'.", 
                                                tags$b("b, "), 
                                                "Plot displaying how statistical power 
                                           varies with maximum follow-up time. The dashed line indicates the statistical power threshold of 0.8.",
                                                tags$b("c, "),
                                                "Side panel for performing power analysis by varying the frequency of measurements. Note that the type of analysis is set to 'freq.'",
                                                tags$b("d, "),
                                                "Plot displaying how statistical power varies with the frequency of
                                           measurements. The dashed line indicates the statistical power threshold of 0.8.")
                                
                              ),
                              tags$p("In Fig. 10b, it can be seen how the power varies with the maximum time of follow-up. As expected, when the follow-up period is short,",
                                     "the statistical power is low, and it increases as the follow-up time increases, until it reaches a point where",
                                     "the power does not increase significantly."),
                              tags$h5("Frequency of Measurements"),
                              tags$p("Another question that researchers may have is how frequently the measurements should be taken. Of course, this would depend, among other factors,",
                                     "on practical considerations, but it can be simulated to estimate the ideal frequency of measurements."),
                              tags$p("Again, for this analysis, all the parameters are fixed, and the only variation is how frequently the measurements are taken,",
                                     "which is defined by the number of evenly spaced measurements performed during the follow-up period.",
                                     "For this example, the analysis evaluated how the power varied from 3 to 26 evenly spaced measurements, with 5 subjects per group and a follow-up",
                                     "period of 42 days (Fig. 10c)."),
                              tags$p("As shown in Fig. 10d, the power increases as the number of measurements increases, although the difference is not dramatic.")
                              ),
                         
                         card(tags$h4(id = "variability", "Data Variability Power Analysis"),
                              tags$p("Finally, another factor that can affect statistical power is the variability in the data. The variability in the model",
                                     "is represented by the standard deviation of the random effects (between-subject variability) and the standard",
                                     "deviation of the residuals (within-subject variability), which can be obtained from the table shown in Fig. 1c."),
                              tags$p("Additionally, the statistical power is also influenced by the magnitude of the differences: the greater the effect",
                                     "of the drug combination, the higher the statistical power to detect synergistic effects."),
                              tags$p("The user can test how statistical power varies by modifying the variability and the drug combination effects in the",
                                     "'Variability Power Analysis' side panel (Fig. 11a)."),
                              tags$p("In this case, the parameters that vary are the standard deviations of the random effects and the residuals, as",
                                     "well as the effect of the drug combination, represented by the coefficient of the drug combination group",
                                     "(i.e., tumor growth rate in the drug combination group). By default, these parameters are evaluated over a range",
                                     "from 10% to 200% of the values obtained from the model fitted in the ", tags$b("Model Estimation"),
                                     "tab, but they can be modified by checking the 'Show Advanced Options' box."),
                              
                              tags$figure(
                                tags$img(
                                  width = "100%",
                                  src = "Fig11.png"
                                ),
                                tags$figcaption(tags$b("Figure 11. a, "), 
                                                "Side panel for performing power analysis by varying the variability and the effect of the drug combination.", 
                                                tags$b("b, "), 
                                                "Plots showing the regression lines and values of the exemplary data, along with the statistical power as a 
                                                function of the standard deviation (SD) of the random effects and residuals, as well as the coefficient for the
                                                drug combination group. The dashed line indicates the statistical power threshold of 0.8.",
                                                "The red dot and vertical dashed line indicate the power for the given SD of random effects and residuals, and growth rate given by the original model.")
                              ),
                              tags$p("From the results shown in Fig. 11b, it can be observed how the power increases as the variability",
                                     "decreases and as the coefficient for the drug combination group decreases (i.e., the growth rate",
                                     "for the drug combination is smaller, indicating a higher drug combination effect).")
                              ),
                         
                    ),
                    card(tags$b("Note: All a priori analyses are conducted with certain parameters held constant while varying others.",
                                "This approach provides users with insights into optimal experimental conditions for achieving adequate",
                                "statistical power. However, this method may not accurately reflect real-world scenarios, where changes",
                                "to one variable can influence others. For instance, increasing or decreasing the sample size or the",
                                "frequency of measurements can impact the variability in the data."))
                    ),
               card(id = "workflow",
                    card_header(h2(
                      "Recommended workflow for using SynergyLMM"
                    )),
                    tags$figure(
                      tags$img(width = "75%", src = "Flowchart_v2.png", alt = "SynergyLMM workflow overview."),
                      tags$figcaption(tags$br(), tags$b("Recommended workflow for using SynergyLMM to assess drug combination effects."),
                                      "HSA: highest single agent. RA: response additivity.")
                      
                    ))
               )
             )
    ),
    
    
    # 1. Tab for the tumor growth analysis ----
    tabPanel("Model Estimation",
             sidebarLayout(
               sidebarPanel(
                 h5("Upload Data") %>% helper(type = "markdown", content = "data_upload"),
                 tags$a(href = "eg_grwth_data.xlsx", download = NA, "Download Example Dataset"),
                 tags$br(),
                 helpText("Supported file formats: TXT, CSV, XLSX."),
                 fileInput("file", "Choose a file (TXT, CSV, XLSX)", 
                           accept = c(".txt", ".csv", ".xlsx")),
                 uiOutput("sample_selector"),   # UI for selecting the Sample column
                 uiOutput("time_selector"),     # UI for selecting the Time column
                 uiOutput("treatment_selector"),# UI for selecting the Treatment column
                 uiOutput("tv_selector"),       # UI for selecting the TV column
                 helpText("Set to minimum value of 'Time' column by default."),
                 numericInput("time_start", "Time Start (First Time Point):", value = NA) %>% helper(type = "markdown", content = "time_start"),
                 helpText("Set to maximum value of 'Time' column by default."),
                 numericInput("time_end", "Time End (Last Time Point):", value = NA) %>% helper(type = "markdown", content = "time_end"),
                 uiOutput("control_group"),
                 uiOutput("drug_a"),
                 uiOutput("drug_b"),
                 uiOutput("drug_c"),
                 uiOutput("combination"),
                 numericInput("min_observations", "Minimum Observations per Subject:", value = 1) %>% helper(type = "markdown", content = "min_obs"),
                 selectInput("grwth_model", label = "Growth Model", choices = list("Exponential" = "exp", "Gompertz" = "gompertz"), selected = "exp") %>%
                   helper(type = "markdown", content = "grwth_model", size = "l"),
                 checkboxInput("show_plot", "Show Plot", value = TRUE),
                 checkboxInput("model_advanced", "Show Advanced Options", value = FALSE),
                 conditionalPanel(
                   condition = "input.model_advanced == true",
                   h5("Define unequal variances for the errors") %>% helper(type = "markdown", content = "unequal_var"),
                   helpText("Check this option to allow the within-group errors have",
                            "unequal variances."),
                   checkboxInput("model_heteroscedasticity", "Describe within-group heteroscedasticity structure.", value = FALSE),
                   conditionalPanel(
                     condition = "input.model_heteroscedasticity == true",
                     helpText("Select 'SampleID', 'Time' or 'Treatment' to use a different variance per subject",
                              "time point or treatment, respectively"),
                     selectInput("varIdent", "Variance Function:", 
                                 choices = c("SampleID", "Time", "Treatment"), selected = "SampleID")),
                   h5("Transformation for null tumor measurements") %>% helper(type = "markdown", content = "tum_vol_0"),
                   helpText("Possible options are 'ignore' (default), to ignore these measurements, or 'transform', to add",
                            "1 unit to all."),
                   selectInput("tum_vol_0", "Null Tumor Measurement Transformation", 
                               choices = c("ignore", "transform"), selected = "ignore"),
                   helpText("Use 'Define' to specify the start values, or 'selfStart' for self-starting values."),
                   conditionalPanel(
                     condition = "input.grwth_model == 'gompertz'",
                     selectInput("selfStart_sel", label = "Start values for Gompertz model", choices = c("selfStart", "Define"), selected = "Define"),
                     conditionalPanel(
                       condition = "input.selfStart_sel == 'Define'",
                       numericInput("r0", "r0 value", value = 0.05),
                       numericInput("rho", "rho value", value = 0.01)
                     )
                   ),
                   helpText("Increasing these values can help when convergence",
                            "errors appear."),
                   h5("Control values for the estimation algorithm") %>% helper(type = "markdown", content = "lmecontrol"),
                   helpText("Maximum number of iterations for the optimization algorithm. Default is 50."),
                   numericInput("maxIter", "maxIter", value = 50),
                   helpText("Maximum number of iterations for the optimization step inside the optimization. Default is 50."),
                   numericInput("msMaxIter", "msMaxIter", value = 50),
                   helpText("Number of iterations for the expectation-maximization algorithm used to refine the initial estimates of the random effects variance-covariance coefficients. Default is 25."),
                   numericInput("niterEM", "niterEM", value = 25),
                   conditionalPanel(
                     condition = "input.grwth_model == 'exp'",
                     helpText("Maximum number of evaluations of the objective function permitted for 'nlminb' optimizer. Default is 200."),
                     numericInput("msMaxEval", "msMaxEval", value = 200)
                   ),
                   conditionalPanel(
                     condition = "input.grwth_model == 'gompertz'",
                     helpText("maximum number of iterations for the PNLS optimization step inside the nlme optimization. Default is 7."),
                     numericInput("pnlsMaxIter", "pnlsMaxIter", value = 7)
                   ),
                   helpText("The optimizer to be used, with 'nlminb' as the default."),
                   conditionalPanel(
                     condition = "input.grwth_model == 'exp'",
                     selectInput("opt", "opt", 
                                 choices = c("nlminb", "optim"), selected = "nlminb")
                   ),
                   conditionalPanel(
                     condition = "input.grwth_model == 'gompertz'",
                     selectInput("opt", "opt", 
                                 choices = c("nlminb", "nlm"), selected = "nlminb")
                   ),
                 ),
                 actionButton("run_analysis", "Run Analysis")
                 ), 
               mainPanel(
                 card(h6("Input Table")),
                 DTOutput("dataTable"),
                 card(h3(textOutput("model_estimates_header")),
                 tableOutput(outputId = "model_estimates")),
                 hidden(downloadButton("downloadEstimates", "Download Model Estimates Table")),
                 hidden(selectInput("estimates_file_type", "Select file type:",
                             choices = c("CSV" = "csv", "TXT" = "txt", "XLSX" = "xlsx"))),
                 card(h3(textOutput("model_plot_header")),
                 plotly::plotlyOutput("model_plot", width = "100%", height = 800)),
                 hidden(downloadButton("downloadModelPlot", "Download Plot")),
                 hidden(selectInput("model_plot_file", "Select file type:",
                                    choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg")))
                 
               )
             )
    ),
    
    # 2. Second tab for synergy calculation ----
    tabPanel("Synergy Analysis",
             sidebarLayout(
               sidebarPanel(
                 h3("Synergy Analysis"),
                 helpText("Click on 'Run Synergy Calculation' to calculate synergy based on the selected reference model."),
                 selectInput("synergy_method", "Select Synergy Calculation Reference Model", 
                             choices = c("Bliss", "HSA", "RA")) %>% helper(type = "markdown", content = "synergy_model"),
                 checkboxInput("synergy_advanced", "Show Advanced Options", value = FALSE),
                 conditionalPanel(
                   condition = "input.synergy_advanced == true",
                   helpText("Select the earliest time point from which to start calculating synergy."),
                   numericInput("min_time", "Minimum Time for Synergy Calculation", value = 0),
                   helpText("Confidence level (1-alpha) to use for confidence intervals and obtain p-values. Must be in the range (0,1)."),
                   numericInput("conf_level", "Confidence level", value = 0.95, min = 0, max = 1),
                   selectInput("padj", "Correction method for adjusting p-values", 
                               choices = c("none", "BH","fdr", "holm", "hochberg", "hommel", "bonferroni", "BY")) %>% helper(type  = "markdown", content = "padj", size = "l"),
                   helpText("Use Robust Estimators to deal with possible model misspecifications."),
                   conditionalPanel(
                     condition = "input.grwth_model == 'exp'",
                     checkboxInput("robustSE", "Use Robust Estimators", value = TRUE) %>% helper(type = "markdown", content = "robust_estimators", size = "l"),
                     selectInput("se_type", "Select which small-sample adjustment should be used", 
                                 choices = c("CR2", "CR0", "CR1", "CR1p", "CR1S", "CR3"))
                   ),
                   conditionalPanel(
                     condition = "input.grwth_model == 'gompertz'",
                     checkboxInput("robustSE", "Use Robust Estimators", value = FALSE) %>% helper(type = "markdown", content = "robust_estimators", size = "l"),
                     selectInput("se_type", "Select which small-sample adjustment should be used", 
                                 choices = c("CR2", "CR0", "CR1", "CR1p", "CR1S", "CR3"))
                   ),
                   helpText("Check this option to fix the seed for the simulation for RA reference model",
                            "and obtain reproducible results."),
                   checkboxInput("setseed_ra", "Obtain reproducible results", value = FALSE),
                   helpText("Set the number of random resampling to calculate the synergy for Response Additivity model."),
                   numericInput("nsim", "Number of simulations for RA model or synergy calculation
                                using the Gompertz growth model", value = 1000)
                 ),
                 actionButton("syn_calc", "Run Synergy Calculation"), width = 3,
               ),
               mainPanel(
                 h3(textOutput("synergy_plot_header")),
                 br(),
                 plotOutput("synergy_plot", height = 800),
                 br(),
                 hidden(downloadButton("downloadSynPlot", "Download Plot")),
                 hidden(selectInput("syn_plot_file", "Select file type:",
                                    choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"))),
                 h3(textOutput("synergy_results_table_header")),
                 br(),
                 DTOutput("synergy_results_table"),
                 hidden(downloadButton("downloadSynergy", "Download Synergy Results")),
                 hidden(selectInput("synergy_file_type", "Select file type:",
                                    choices = c("CSV" = "csv", "TXT" = "txt", "XLSX" = "xlsx"))),
                 h3(textOutput("synergy_estimates_table_header")),
                 br(),
                 DTOutput("synergy_estimates_table"),
                 hidden(downloadButton("downloadSynEstimates", "Download Model Estimates")),
                 hidden(selectInput("syn_estimates_file_type", "Select file type:",
                                    choices = c("CSV" = "csv", "TXT" = "txt", "XLSX" = "xlsx")))
               )
             )
    ),
    
    # 3. Third tab for model diagnostics ----
    tabPanel("Model Diagnostics",
             sidebarLayout(
               sidebarPanel(
                 card(h3("Model Diagnostics") %>% helper(type = "markdown", content = "model_diag", size = "l"),
                 helpText("Plots and results shown in this tab help verify that the model is appropriate for the data.",
                          "Reliable results from statistical tests require that model assumptions — such as normality and homoscedasticity of random effects and residuals — are met.",
                          "Therefore, users are recommended to thoroughly evaluate the estimated model, as certain experimental datasets may result in highly imprecise estimates.",
                          "More information about interpretation of model diagnostics can be found clicking in the '?' icon.")),
                 card(h3("Influential Diagnostics") %>% helper(type = "markdown", content = "influential_diag", size = "l"),
                 helpText("Influential diagnostics allow the user to identify highly influential subjects (animals) in the experiment.", 
                          "It is based on detecting those subjects with a higher influence on the fitted values, or in the estimation of control and",
                          "treatment groups estimated coefficients. In both cases, the results are evaluated using Cook’s distances.",
                          "More information about interpretation of model diagnostics can be found clicking in the '?' icon.")),
                 checkboxInput("influential_advanced", "Show Advanced Options", value = FALSE),
                 conditionalPanel(
                   condition = "input.influential_advanced == true",
                   selectInput("diag_norm_test", 
                               label = "Test for normality assessment:", 
                               choices = c("Shapiro - Wilk Normality Test" = "shapiroTest",
                                           "D'Agostino Normality Test" = "dagoTest",
                                           "Anderson - Darling Normality Test" = "adTest"),
                               selected = "shapiroTest"),
                   helpText("Threshold for the p-value of outlier observations based on their normalized residuals."),
                   numericInput("diag_pval",
                                label = "Threshold p-value for outliers:",
                                value = 0.05, min = 0, max = 1),
                   helpText("Set a numeric value indicating the threshold for the Cook's distance. If not specified, the threshold is set to three times the mean of the Cook's distances values."),
                   numericInput("cook_thr", "Threshold for Cook's distance:", value = NA),
                   helpText("Type of Cook's distance to calculated. Possible options are 'fitted', to calculte Cook's distances based on the change in fitted values, 
                            or 'fixef' to calculate Cook's distances based on the change in the fixed effects. "),
                   selectInput("cook_type", label = "Type of Cook's distance to calculated:",
                               choices = c("fitted", "fixef"), selected = "fitted"),
                   helpText("Cook's distances calculation implies re-fitting multiple models. In some cases,
                            the fitting process can fail due to convergence issues, especially when using the Gompertz model. 
                            The maximum number of iterations can be increased to solve those issues."),
                   numericInput("maxIter_cook", label = "Limit of maximum number of iterations for the optimization algorithm:",
                                value = 1000)
                 ),
                 width = 3,
               ),
               mainPanel(
                 card(h1(textOutput("model_diagnostics_header")),
                 plotOutput("model_diag_plot", height = 1200)),
                 hidden(downloadButton("downloadDiagPlot", "Download Plot")),
                 hidden(selectInput("diag_plot_file", "Select file type:",
                                    choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"))),
                 h3(textOutput("normality_ranef_header")),
                 verbatimTextOutput("ranef_norm"),
                 h3(textOutput("normality_resid_header")),
                 verbatimTextOutput("resid_norm"),
                 card(h3(textOutput("outliers_header")), DTOutput("resid_outliers")), 
                 hidden(downloadButton("downloadOutliers", "Download Outliers Table")),
                 hidden(selectInput("outliers_file_type", "Select file type:",
                                    choices = c("CSV" = "csv", "TXT" = "txt", "XLSX" = "xlsx"))),
                 h3(textOutput("cook_dist_header")),
                 plotOutput("cook_dist", height = 800),
                 verbatimTextOutput("cook_dist_sub"),
                 hidden(downloadButton("downloadCooksDist", "Download Cook's Distances")),
                 br()
               )
             )
    ),
    
    # 4. Fourth tab for model performance ----
    tabPanel("Model Performance",
             sidebarLayout(
               sidebarPanel(
                 h3("Model Performance") %>% helper(type = "markdown", content = "ObsvsPred"),
                 numericInput("nrow", "Number of rows for ploting Observed vs Predicted results:", value = 4),
                 numericInput("ncol", "Number of columns for ploting Observed vs Predicted results:", value = 5),
                 width = 3,
               ),
               mainPanel(
                 card(h3(textOutput("performance_header")),
                 tableOutput(outputId = "model_performance")),
                 h3(textOutput("obs_vs_plot_header")),
                 plotOutput("obs_vs_plot", height = 1200),
                 hidden(downloadButton("downloadObsPredPlot", "Download Plot")),
                 hidden(selectInput("obs_pred_file", "Select file type:",
                                    choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg")))
               )
             )
    ),
    
    # 5. Fifth tab fo power calculation ----
    tabPanel("Power Analysis",
             sidebarLayout(
               sidebarPanel(
                 h2("Post Hoc Power Analysis") %>% helper(type = "markdown", content = "posthoc"),
                 helpText("Click on 'Run Post Hoc Power Analysis' to calculate the statistical power for the current dataset."),
                 selectInput("pwr_method", "Select Synergy Calculation Reference Model", 
                             choices = c("Bliss", "HSA")),
                 uiOutput("posthocpwr_time"),
                 checkboxInput("posthocpwr_advanced", "Show Advanced Options", value = FALSE),
                 conditionalPanel(condition = "input.posthocpwr_advanced == true",
                                  numericInput("pwr_nsim", "Number of simulations to run:", value = 1000),
                                  helpText("Check this option to fix the seed for the simulation",
                                           "and obtain reproducible results."),
                                  checkboxInput("setseed", "Obtain reproducible results.", value = FALSE),
                                  numericInput("pwr_pval", "p-value threshold for power calculation:", value = 0.05)
                 ),
                 actionButton("post_hoc_pwr", "Run Post Hoc Power Analysis"),
                 tags$hr(), # horizontal line
                 tags$br(), # line break
                 
                 tags$h2("A Priori Power Analysis") %>% helper(type = "markdown", content = "apriori", size = "l"),
                 tags$h4("Sample Size Power Analysis"),
                 helpText("Click on 'Run Sample Size Power Analysis' to calculate the statistical 
                          power for a model with the calculated estimates from the actual dataset varying the sample size."),
                 sliderInput("sample_n", "Sample size per group:", min = 0, max = 100, step = 1, value = c(5,10)),
                 selectInput("sample_size_pwr_method", "Select Synergy Calculation Reference Model", 
                             choices = c("Bliss", "HSA")),
                 actionButton("sample_size_pwr", "Run Sample Size Power Analysis"), width = 3,
                 tags$hr(), # horizontal line
                 tags$br(), # line break
                 tags$h4("Time Power Analysis"),
                 helpText("Click on 'Run Time Power Analysis' to calculate the statistical 
                          power for a model with the calculated estimates from the actual dataset varying the time of follow-up or the frequency of the measurements."),
                 # Dropdown to choose between "max" or "freq"
                 selectInput("max_freq", "Select analysis based on maximum time of follow-up ('max') or frequency of measurements ('freq')", 
                             choices = c("max", "freq")),
                 
                 # Conditional Panel for "max" selection
                 conditionalPanel(
                   condition = "input.max_freq == 'max'",
                   numericInput("time_pwr_npg", "Sample size per group:", value = 5),
                   numericInput("time_pwr_min_t", "Minimum time of follow-up:", value = 3),
                   numericInput("time_pwr_max_t", "Maximum time of follow-up:", value = 15),
                   numericInput("time_pwr_step", "Time interval between measurements:", value = 3),
                   selectInput("time_pwr_method", "Select Synergy Calculation Reference Model", 
                               choices = c("Bliss", "HSA")),
                 ),
                 # Conditional Panel for "freq" selection
                 conditionalPanel(
                   condition = "input.max_freq == 'freq'",
                   numericInput("time_pwr_npg", "Sample size per group:", value = 5),
                   numericInput("time_pwr_min_freq", "Minimum number of measurements:", value = 3),
                   numericInput("time_pwr_max_freq", "Maximum number of measurements:", value = 5),
                   numericInput("time_pwr_max_t", "Time of follow-up:", value = 15),
                   selectInput("time_pwr_method", "Select Synergy Calculation Reference Model", 
                               choices = c("Bliss", "HSA")),
                 ),
                 actionButton("time_pwr", "Run Time Power Analysis"),
                 tags$hr(), # horizontal line
                 tags$br(), # line break
                 tags$h4("Variability Power Analysis"),
                 helpText("Click on 'Run Variability Power Analysis' to calculate the statistical", 
                          "power for a model with especific values estimates with different",
                          "variability and/or different effect size of the drug combination.",
                          "(If the values for the estimates are not especified, the values from the model are used by default)."),
                 numericInput("var_pwr_npg", "Sample size per group:", value = 5),
                 selectInput("var_pwr_method", "Select Synergy Calculation Reference Model", 
                             choices = c("Bliss", "HSA")),
                 checkboxInput("var_pwr_check", "Show Advanced Options", value = FALSE),
                 helpText("Use advanced options to set specific values of the estimates",
                          "to test the power, such us changing the growth rate for the",
                          "drug combination", "or increasing/decreasing the variance of the residuals."),
                 conditionalPanel(
                   condition = "input.var_pwr_check == true",
                   h5("Range of variability and/or drug combination effect size to calculate power."),
                   helpText("If not especified, the minimum and maximum values are",
                            "set to 10% and 200% of the model estimates by default."),
                   numericInput("sd_eval_min", "Minimum value for random effect standard deviation to test:", value = NA),
                   numericInput("sd_eval_max", "Maximum value for random effect standard deviation to test:", value = NA),
                   numericInput("sgma_min", "Minimum value for residuals standard deviation to test:", value = NA),
                   numericInput("sgma_max", "Maximum value for residuals standard deviation to test:", value = NA),
                   numericInput("grwrComb_min", "Minimum value for growth rate for Drug A+B Combination group:", value = NA),
                   numericInput("grwrComb_max", "Maximum value for growth rate for Drug A+B Combination group:", value = NA),
                   h5("Estimated parameters from the model"),
                   helpText("If not specified, the values for the current model are used by default."),
                   numericInput("grwrControl_pwr_var", "Estimated growth rate for Control group:", value = NA),
                   numericInput("grwrA_pwr_var", "Estimated growth rate for Drug A group:", value = NA),
                   numericInput("grwrB_pwr_var", "Estimated growth rate for Drug B group:", value = NA),
                   numericInput("grwrComb_pwr_var", "Estimated growth rate for Drug A+B Combination group:", value = NA),
                   numericInput("sd_ranef", "Random effects standard deviation (between-subject variance) for the model:", value = NA),
                   numericInput("sgma", "Residuals standard deviation (within-subject variance) for the model:", value = NA)
                   ),
                 tags$br(),
                 actionButton("var_pwr", "Run Variability Power Analysis")
               ),
               mainPanel(
                 br(),
                 h3(textOutput("post_hoc_pwr_header")),
                 card(textOutput("post_hoc_pwr")),
                 br(),
                 card(
                 h3(textOutput("sample_size_pwr_header")),
                 plotOutput("sample_size_pwr")),
                 hidden(downloadButton("downloadSampleSizePwr", "Download Plot")),
                 hidden(selectInput("size_pwr_plot_file", "Select file type:",
                                    choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"))),
                 br(),
                 card(h3(textOutput("time_pwr_header")),
                 plotOutput("time_pwr")),
                 hidden(downloadButton("downloadTimePwr", "Download Plot")),
                 hidden(selectInput("time_pwr_plot_file", "Select file type:",
                                    choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg"))),
                 br(),
                 card(h3(textOutput("var_pwr_header")),
                 plotOutput("var_pwr")),
                 hidden(downloadButton("downloadVarPwr", "Download Plot")),
                 hidden(selectInput("var_pwr_plot_file", "Select file type:",
                                    choices = c("PDF" = "pdf", "PNG" = "png", "SVG" = "svg")))
               )
             )
    ),
     # 6. Citation tab ----
    tabPanel(
      title = "Citation",
      fluidPage(
        tags$br(),
        h2("How to cite this work"),
        h5(
          "If you use this app or the underlying R package in your work, ",
          "please cite the following article:"
        ),
        
        ## Primary journal article
        card(
        h3("Primary article"),
        h5("Romero-Becerra, R., Zhao, Z., Nebdal, D.", tags$i("et al."), 
               tags$b("Improved analysis of in vivo drug combination experiments with a comprehensive statistical framework and web-tool."),
               tags$i("Nat Commun"),tags$b("16"),", 10151 (2025).",
               tags$a(href = "https://doi.org/10.1038/s41467-025-65218-9", "https://doi.org/10.1038/s41467-025-65218-9")),
        ),
        
        ## Links
       card(
         h3("Other Resources"),
         tags$ul(
           tags$li(
             tags$a(href="https://cran.r-project.org/package=SynergyLMM", "CRAN package page")
           ),
           tags$li(
             tags$a(href = "https://github.com/RafRomB/SynergyLMM", "GitHub repository")
           ),
           tags$li(
             tags$a(href = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-65218-9/MediaObjects/41467_2025_65218_MOESM1_ESM.pdf#page=1.10",
                    "Guidelines for the correct use of the framework")
           )
         )
       )
      )
    )
  )
)