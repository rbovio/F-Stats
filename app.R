library(shiny)
library(ggplot2)

ui <- fluidPage(
  #sidebarLayout
  sidebarLayout(
    sidebarPanel(style = "position:fixed;width:30%;",
      fluidRow(align = 'center',
               helpText("Choose the number of populations (min=2, max=6) and the size of each population (min=10, max=10,000) for anaylsis"),
               numericInput("num_pops", "number of populations", value = 3, min = 2, max = 6, width = '30%')
      ),
      fluidRow(align = 'center',
               htmlOutput("pop_slider_setup")
      )
    ),#end of sidebarPanel
    
    mainPanel(width = 8,
      h2("F-statistics made EASY!"),
      #create info section
      fluidRow(align = 'justify',
               htmlOutput("fstats_notes"),
      ), 
      #observed genotype counts table
      h3("observed genotype counts", align = 'center'),
      fluidRow(align = 'center',
               tableOutput("observed_genotype_counts_table"),
               plotOutput("observed_genotype_counts_plot")
      ),
      #Step 1: Calculate the allele frequencies for each population
      h3("Step 1: Calculate the allele frequencies for each population", align = 'center'),
      fluidRow(column(width = 6, align = 'right',
                      tableOutput("allele_freq_table")
      ),
      column(width = 6, align = 'left',
             htmlOutput("allele_freq_formula"),
             tags$head(tags$style("#allele_freq_formula{font-size: 15px;
                                 }"
             ))
      )
      ),#end of step 1
      
      #Step 2: Calculate expected genotype counts under Hardy-Weinberg Equilibrium
      h3("Step 2: Calculate expected genotype counts under Hardy-Weinberg Equilibrium", align = 'center'),
      fluidRow(column(width = 6, align = 'right',
                      tableOutput("expected_genotype_counts_table")
      ),
      column(width = 6, align = 'left',
             htmlOutput("expected_genotype_counts_formula"),
             tags$head(tags$style("#expected_genotype_counts_formula{font-size: 15px;
                                   }"
             ))
      )
      ),#end of step 2
      
      #Step 3: Calculate observed and expected heterozygote frequencies and inbreeding coefficients
      h3("Step 3: Calculate observed and expected heterozygote frequencies and local inbreeding coefficients", align = 'center'),
      htmlOutput("step3_notes", align = 'center'),
      tags$head(tags$style("#step3_notes{font-size: 15px;
                                         }"
      )),
      fluidRow(column(width = 6, align = 'right',
                      tableOutput("het_freq_table")
      ),
      column(width = 6, align = 'left',
             htmlOutput("het_freq_formula"),
             tags$head(tags$style("#het_freq_formula{font-size: 15px;
                                         }"
             )),
      )
      ),
      #h4("Step 3 notes: Positive F means fewer heterozygotes than expected; indicates inbreeding. Negative F means more heterozygotes than expected; excess outbreeding. Expected heterozygosity = 0 when one of the alleles is fixed.", align = 'center'),
      #end of step 3
      
      #Step 4: Calculate the global allele frequency
      h3("Step 4: Calculate the global allele frequency", align = 'center'),
      fluidRow(column(width = 6, align = 'right',
                      tableOutput("global_allele_freq_table")
                      ),
               column(width = 6, align = 'left',
                      htmlOutput("global_allele_freq_formula1"),
                      tags$head(tags$style("#global_allele_freq_formula1{font-size: 15px;
                                         }"
                      )),
                      htmlOutput("global_allele_freq_formula2"),
                      tags$head(tags$style("#global_allele_freq_formula2{font-size: 15px;
                                         }"
                      )),
                      )
      ),#end of step 4
      
      #Step 5: Calculate the global heterozygote indicies
      h3("Step 5: Calculate the global heterozygote indicies", align = 'center'),
      fluidRow(column(width = 6, align = 'right',
                      
                      tableOutput("global_het_indicies_table")
                      ),
               column(width = 6, align = 'left',
                      htmlOutput("global_het_indicies_formula_HI_1"),
                      tags$head(tags$style("#global_het_indicies_formula_HI_1{font-size: 15px;
                                         }"
                      )),
                      htmlOutput("global_het_indicies_formula_HI_2"),
                      tags$head(tags$style("#global_het_indicies_formula_HI_2{font-size: 15px;
                                         }"
                      )),
                      htmlOutput("global_het_indicies_formula_HS_1"),
                      tags$head(tags$style("#global_het_indicies_formula_HS_1{font-size: 15px;
                                         }"
                      )),
                      htmlOutput("global_het_indicies_formula_HS_2"),
                      tags$head(tags$style("#global_het_indicies_formula_HS_2{font-size: 15px;
                                         }"
                      )),
                      htmlOutput("global_het_indicies_formula_HT_1"),
                      tags$head(tags$style("#global_het_indicies_formula_HT_1{font-size: 15px;
                                         }"
                      )),
                      htmlOutput("global_het_indicies_formula_HT_2"),
                      tags$head(tags$style("#global_het_indicies_formula_HT_2{font-size: 15px;
                                         }"
                      ))
                      )
               
      ),
      #h4("Step 4 & 5 notes: You need to multiple by the size of the population (N) to account for differences in population sizes. Hi is based on observed heterozygosities in individuals in subpopulations. Hs is based on expected heterozygosities in subpopulations. Ht is based on expected heterozygosities for overall total population", align = 'center'),
      #end of step 5
      
      #Step 6: Calculate the global fixation indicies
      h3("Step 6: Calculate the global fixation indicies", align = 'center'),
      fluidRow(column(width = 6, align = 'right',
                      tableOutput("global_fstats_table")
                      ),
      column(width = 6, align = 'left',
             htmlOutput("global_fstats_formula"),
             tags$head(tags$style("#global_fstats_formula{font-size: 15px;
                                         }"
             ))
             )
      ),
      #h4("Step 6 notes: Fis (inbreeding coefficient of an individual (I) relative to the subpopulation (S)), Fst (inbreeding coefficient of subpopulations (S) relative to the total population (T)), and Fit (inbreeding coefficient of an individual (I) relative to the total (T) population) are related to the amounts of heterozygosity at various levels of population structure. Together, they are called F-statistics, and are derived from F, the inbreeding coefficient. Compare and contrast the global Fis above with the local inbreeding coefficient (F) of Step 3. Here, we are using a weighted average of the individual heterozygosities over all the subpopulations. Both FIS and  Fs are based on the observed heterozygosities, whereas Fst and Fit are based on expected heterozygosities.", align = 'center')
      #end of step 6
      
    )#end of mainPanel
    
  )#end of sidebarLayout

)#end of UI


# Define server logic required
server <- function(input, output) {
  
  #fstat_notes
  output$fstats_notes = renderUI({
    str1 = "One of the results of population structure is a reduction in heterozygosity. When populations split, alleles have a higher chance of reaching fixation within subpopulations, especially if the subpopulations are small or have been isolated for long periods. This reduction in heterozygosity can be thought of as an extension of inbreeding, with individuals in subpopulations being more likely to share a recent common ancestor. F-Statistics use inbreeding coefficients to describe the partitioning of genetic variation within and among populations and can be calculated at three different levels. The first F-statistic, FIS, measures the degree of inbreeding within individuals relative to the rest of their subpopulation. This reflects the probability that two alleles within the same individual are identical by descent, and is the same as the inbreeding coefficient F. It is calculated as:"
    str4 = "FIS = HS - HI / HS"
    str5 = "HI is the observed heterozygosity in a subpopulation at the time of investigation (individual heterozygosity) and HS is the heterozygosity that would be expected if the subpopulation was in HWE (subpopulation heterozygosity)."
    str6 = "The second F-statistic is FST (also known as the fixation index), and this provides an estimate of the genetic differentiation between subpopulations. It is a measure of the degree of inbreeding within a subpopulation relative to the total population (total population here meaning all of the subpopulations combined), and reflects the probability that two alleles drawn at random from within a subpopulation are identical by descent. It is calculated as:"
    str7 = "FST = HT - HS / HT"
    str8 = "where HS is the same as above and HT is the expected heterozygosity of the total population."
    str9 = "The third F-statistic, which is used much less frequently than the other two, is FIT. This provides an overall inbreeding coefficient for an individual by measuring the heterozygosity of an individual relative to the total population. FIT is therefore influenced by both non-random mating within a subpopulation (FIS) and population subdivision (FST), and is calculated as:"
    str10 = "FIT = HT - HI / HT"
    str11 = "where HT and HI are the same as above. The relationship between the three statistics is given as:"
    str12 = "FIT = FIS + FST — (FIS)(FST)"
    str13 = "Since FST measures the extent to which populations have differentiated from one another, this is the F-statistic with which we are most concerned."
    brk = "<br/>"
    HTML(paste(str1,brk,str4,brk,str5,str6,brk,str7,brk,str8,str9,brk,str10,brk,str11,brk,str12,brk,str13, sep = "<br/>"))
  })

  #create output panels for each population
  output$pop_slider_setup = renderUI({
    col_list <- tagList()
    #create a list to store inputs for population size (numericInput) and initial genotype values (uiOutput).
    for(i in 1:input$num_pops){
      col_list[[i]] = column(width = 4, align = 'center',
                             numericInput(paste0("max",i), paste0("pop",i," size"), value = 100, min = 10, max = 10000, width = '100%'),
                             wellPanel(
                               uiOutput(paste0("slider1",i)),
                               uiOutput(paste0("slider2",i))
                             ),
                             tableOutput(paste0("restable",i))
      )
    }
    #Display objects in list to UI
    tagList(
      col_list
    )
  })#end of pop_slider_setup
  
  for(i in 1:10){
    #slider1
    eval(parse(text = (paste0("output$slider1",i," <- renderUI ({sliderInput(\"slider1",i,"\", \"AA\", min = 0,  max = input$max",i,", value = 0, step = 1)})"))))
    #slider2
    eval(parse(text = (paste0("output$slider2",i," <- renderUI ({sliderInput(\"slider2",i,"\", \"Aa\", min = 0,  max = input$max",i," - input$slider1",i,", value = 0, step = 1)})"))))
  }#end of slider for loop
  
  #Step 0. Display observed genotype counts table
  output$observed_genotype_counts_table = renderTable({
    m = matrix(nrow = 4, ncol = input$num_pops)
    str_colnames = c()
    myvals_for_plot = c()
    genos_for_plot = c()
    pops_for_plot = c()
    for(i in 1:input$num_pops){
      eval(parse(text = (paste0("myvals",i," = c(input$slider1",i,", input$slider2",i,", input$max",i,"-input$slider1",i,"-input$slider2",i,", input$max",i,")"))))
      eval(parse(text = (paste0("m[,",i,"] = myvals",i))))
      str_colnames = c(str_colnames,paste0("pop",i))
      myvals_for_plot = c(myvals_for_plot, paste0("myvals",i))
      genos_for_plot = c(genos_for_plot, c("AA", "Aa", "aa"))
      pops_for_plot = c(pops_for_plot, rep(paste0("pop",i),3))
    }
    #table dataframe
    colnames(m) = str_colnames
    m = as.data.frame(m)
    d = data.frame(genotypes = c("AA", "Aa", "aa", "Total"))
    df = cbind(d,m)
    df
  })#end of observed genotype counts table
  
  output$observed_genotype_counts_plot = renderPlot({
    myvals_for_plot = c()
    genos_for_plot = c()
    pops_for_plot = c()
    for(i in 1:input$num_pops){
      eval(parse(text = (paste0("myvals",i," = c(input$slider1",i,", input$slider2",i,", input$max",i,"-input$slider1",i,"-input$slider2",i,")"))))
      myvals_for_plot = c(myvals_for_plot, eval(parse(text = (paste0("myvals",i)))))
      genos_for_plot = c(genos_for_plot, c("AA", "Aa", "aa"))
      pops_for_plot = c(pops_for_plot, rep(paste0("pop",i),3))
    }
    df_plot = data.frame(pop = pops_for_plot,
                         geno = genos_for_plot,
                         values = myvals_for_plot
    )
    #plot dataframe
    p <- ggplot(data = df_plot, aes(x = geno, y = values))
    p <- p + geom_bar(stat = "identity", width = 0.5, position = "dodge")
    p <- p + facet_grid(. ~ pop)
    p <- p + theme_bw()
    p <- p + theme(axis.text.x = element_text(angle = 90))
    p
  })#end of observed genotype counts plot
  
  #Step 1.  Calculate the allele frequencies for each population:
  output$allele_freq_formula = renderUI({
    str1 = paste0("<b><i>p</i><sub>1</sub></b> = (2 x Obs<sub>AA</sub> + Obs<sub>Aa</sub>) / (2 x N)")
    str2 = paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',(2*input$slider11+input$slider21)/(2*input$max1)),"</b>"," = (2 x ",input$slider11," + ",input$slider21,") / (2 x ", input$max1,"</font>)" )
    str3 = "</br>"
    HTML(paste(str1, str2, str3, sep = '</br>'))
  })

  output$allele_freq_table = renderTable({
    m = matrix(nrow = 2, ncol = input$num_pops)
    pvals = c()
    str_colnames = c()
    for(i in 1:input$num_pops){
      pvals = c(pvals,  eval(parse(text = (paste0("(2*input$slider1",i," + input$slider2",i,")/(2*input$max",i,")")))))
      str_colnames = c(str_colnames,paste0("pop",i))
    }
    m[1,] = pvals
    m[2,] = 1-pvals
    colnames(m) = str_colnames
    m = as.data.frame(m)
    d = data.frame(allele = c("p", "q"))
    df = cbind(d,m)
    df
  })#end of allele freq table
  
  #Step 2. Calculate the expected genotypic counts under Hardy-Weinberg Equilibrium, and then calculate the excess or deficiency of homozygotes in each subpopulation.
  #expected_AA_genotype_counts = N*p^2
  #expected_Aa_genotype_counts = N*2*p*q
  #expected_aa_genotype_counts = N*q^2
  #N = input$max
  #p = (2*input$slider1+input$slider2)/(2*input$max)
  #q = 1 - p
  
  output$expected_genotype_counts_formula = renderUI({
    p = (2*input$slider11+input$slider21)/(2*input$max1)
    q = 1-p
    str1 = paste0("<b>Exp<sub>AA,1</sub></b> = N x <i>p</i><sup>2</sup>")
    str2 = paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',input$max1*p^2),"</b>  = ",input$max1," x ", sprintf(fmt = '%#.3f',p), "<sup>","2","</sup></font>")
  
    str3 = paste0("<b>Exp<sub>Aa,1</sub></b> = N x 2<i>p</i><i>q</i>")
    str4 = paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',input$max1*2*p*q),"</b>"," = ",input$max1," x 2(",sprintf(fmt = '%#.3f',p),")(",sprintf(fmt = '%#.3f',q),")","</font>")
    
    str5 = paste0("<b>Exp<sub>aa,1</sub></b> = N x <i>q</i><sup>2</sup>")
    str6 = paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',input$max1*q^2),"</b>  = ",input$max1," x ", sprintf(fmt = '%#.3f',q), "<sup>","2","</sup></font>")
    
    HTML(paste(str1, str2, str3, str4, str5, str6, sep = '</br>'))
  })
  
  output$expected_genotype_counts_table = renderTable({
    m = matrix(nrow = 3, ncol = input$num_pops)
    myvals = c()
    str_colnames = c()
    for(i in 1:input$num_pops){
      eval(parse(text = (paste0("myvals",i," = c( (input$max",i," * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,"))^2)  ,  (input$max",i," * 2 * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")) * (1 - ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")))) ,  (input$max",i," * (1-((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")))^2))"))))
      eval(parse(text = (paste0("m[,",i,"] = myvals",i))))
      str_colnames = c(str_colnames,paste0("pop",i))
    }
    colnames(m) = str_colnames
    m = as.data.frame(m)
    d = data.frame(genotypes = c("AA", "Aa", "aa"))
    df = cbind(d,m)
    df
  })#end of expected genotype counts table
  
  #Step 3. Calculate the local observed and expected heterozygosity of each subpopulation
  #observed_het_freq = input$slider2/input$max
  #for two locus model: expected_het_freq = 2*p*q
  #Step 5. Calculate the local inbreeding coefficient of each subpopulation 
  #positive F means fewer heterozygotes than expected indicates inbreeding
  #negative F means more heterozygotes than expected means excess outbreeding
  #het_exp = 0 when one of the alleles is fixed. 
  #F = (het_exp - het_obs) / het_exp
  
  output$step3_notes = renderUI({
    str1 = paste0("When H<sub>obs,i</sub> = 0 & p > 0 & q > 0 then F = 1")
    str2 = paste0("When H<sub>obs,i</sub> = N<sub>i</sub> then F = -1")
    str3 = paste0("When H<sub>obs,i</sub> = H<sub>exp,i</sub> then F = 0")
    str4 = paste0("positive F means fewer heterozygotes than expected; indicates inbreeding")
    str5 = paste0("negative F means more heterozygotes than expected; indicates outbreeding")
    str6 = paste0("F = 0 means the population is consistent with HWE")
    str7 = "</br>"
    HTML(paste(str1,str4,str2,str5,str3,str6,str7,sep = "</br>"))
  })
  
  output$het_freq_formula = renderUI({
    p = (2*input$slider11+input$slider21)/(2*input$max1)
    q = 1-p
    hobs = input$slider21/input$max1
    hexp = 2*p*q
    str1 = paste0("<b><i>H</i><sub>obs,1</sub></b> = Obs<sub>Aa</sub> / N")
    str2 = paste0("<font color=\"#FF0000\"><b>", sprintf(fmt = '%#.3f',input$slider21/input$max1),"</b> = ", input$slider21," / ", input$max1,"</font>")
    
    str3 = paste0("<b><i>H</i><sub>exp,1</sub></b> = 2<i>pq</i>")
    str4 = paste0("<font color=\"#FF0000\"><b>", sprintf(fmt = '%#.3f',2*p*q),"</b> = 2(", sprintf(fmt = '%#.3f',p), ")(",sprintf(fmt = '%#.3f',q),")","</font>")
    
    str5 = paste0("<b><i>F</i><sub>1</sub></b> = <i>H</i><sub>exp,1</sub> - <i>H</i><sub>obs,1</sub> / <i>H</i><sub>exp,1</sub>")
    str6 = paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',(hexp - hobs) / hexp),"</b> = ", sprintf(fmt = '%#.3f',hexp)," - ", sprintf(fmt = '%#.3f',hobs)," / ",sprintf(fmt = '%#.3f',hexp), "</font>")
    HTML(paste(str1, str2, str3, str4, str5, str6, sep = '</br>'))
    
  })
  
  output$het_freq_table = renderTable({
    m = matrix(nrow = 3, ncol = input$num_pops)
    myvals_obs = c()
    myvals_exp = c()
    str_colnames = c()
    for(i in 1:input$num_pops){
      myvals_obs = c(myvals_obs,  eval(parse(text = (paste0("input$slider2",i," /input$max",i)))))
      myvals_exp = c(myvals_exp,  eval(parse(text = (paste0("2 * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")) * (1 - ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")))")))))
      str_colnames = c(str_colnames,paste0("pop",i))
    }
    m[1,] = myvals_obs
    m[2,] = myvals_exp
    m[3,] = (myvals_exp - myvals_obs) / myvals_exp
    colnames(m) = str_colnames
    m = as.data.frame(m)
    d = data.frame(heterozgosity = c("observed", "expected", "F"))
    df = cbind(d,m)
    df
  })#end of heterozygote frequency table
  
  #Step 4. Calculate (p-bar, the frequency of allele A) and q-bar over the total population.
  #p-bar = (pi*Ni) + (pi+1*Ni+1) / (Ni + Ni+1)
  #Note that we weight by population size. 
  
  output$global_allele_freq_formula1 = renderUI({
    num_formula1 = c()
    den_formula1 = c()
    for(i in 1:input$num_pops){
      num_formula1 = c(num_formula1, paste0("N","<sub>",i,"</sub>"," x p","<sub>",i,"</sub>"), " + ")
      den_formula1 = c(den_formula1, paste0("N","<sub>",i,"</sub>"), " + ")
    }
    #remove the plus sign at end of num_formula1 vector
    num_formula1 = num_formula1[-length(num_formula1)]
    den_formula1 = den_formula1[-length(den_formula1)]
    formula1 = c(paste0("<b>","pbar","</b>"," = "),num_formula1, "/", den_formula1)
    #construct string outputs
    str1 = formula1
    HTML(paste(str1))
  })

  output$global_allele_freq_formula2 = renderUI({
    num_formula2.1 = c()
    num_formula2.2 = c()
    den_formula2 = c()
    for(i in 1:input$num_pops){
      num_formula2.1 = c(num_formula2.1, eval(parse(text = (paste0("input$max",i)))))
      num_formula2.2 = c(num_formula2.2, eval(parse(text = (paste0("((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,"))")))))
      den_formula2 = c(den_formula2, eval(parse(text = (paste0("input$max",i)))))
    }
    #calculate pbar
    pbar = sum(num_formula2.1*num_formula2.2) / sum(den_formula2)
    #add "x" and "+" to num_formula2.1 and num_formula2.2, respectively
    x_var = c(rep(" x ", input$num_pops))
    plus_var = c(rep(" + ", input$num_pops))
    #combine numerator values and operators
    num = c(rbind(num_formula2.1, x_var, sprintf(fmt = '%#.3f',num_formula2.2), plus_var))
    den = c(rbind(den_formula2, plus_var))
    #remove the plus sign at end of num vector
    num = num[-length(num)]
    den = den[-length(den)]
    #construct formula to print
    formula2 = c(paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',pbar),"</font></b>"," = "),"<font color=\"#FF0000\">",num,"</font>","/","<font color=\"#FF0000\">",den,"</font>")
    #construct string outputs
    str1 = formula2
    HTML(paste(str1))

  })
  
  output$global_allele_freq_table = renderTable({
    m = matrix(nrow = 2, ncol = 1)
    num = c()
    den = c()
    for(i in 1:input$num_pops){
      num = c(num, eval(parse(text = (paste0("input$max",i," * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,"))")))))
      den = c(den, eval(parse(text = (paste0("input$max",i)))))
    }
    pbar = sum(num) / sum(den)
    m[1,] = pbar
    m[2,] = 1-pbar
    colnames(m) = "frequency"
    m = as.data.frame(m)
    d = data.frame(allele = c("pbar", "qbar"))
    df = cbind(d,m)
    df
  })#end of global allele freq table
  
  #Step 5. Calculate the global heterozygosity indices (over Individuals, Subpopulations and Total population)
  #need to multiple by population size (N) to account for differences in population sizes 
  #Note that the first two calculations employ a weighted average of the values in the whole set of subpopulations.
  #Hi based on observed heterozygosities in individuals in subpopulations 
  #Hi = ((Hobs1*N1) + (Hobs2*N2)) / Ntotal
  #Hs based on expected heterozygosities in subpopulations
  #Hs = ((Hexp1*N1) + (Hexp2*N2)) / Ntotal
  #Ht based on expected heterozygosities for overall total population
  #Ht = 2*pbar*(1-pbar)
  
  #HI formula 1
  output$global_het_indicies_formula_HI_1 = renderUI({
      num_formula1 = c()
      den_formula1 = c()
      for(i in 1:input$num_pops){
        num_formula1 = c(num_formula1, paste0("N","<sub>",i,"</sub>"," x <i>H</i>","<sub>","Obs",i,"</sub>"), " + ")
        den_formula1 = c(den_formula1, paste0("N","<sub>",i,"</sub>"), " + ")
      }
      #remove the plus sign at end of num_formula1 vector
      num_formula1 = num_formula1[-length(num_formula1)]
      den_formula1 = den_formula1[-length(den_formula1)]
      formula1 = c(paste0("<b><i>H</i><sub>I</sub></b>"," = "),num_formula1, "/", den_formula1)
      #construct string outputs
      str1 = formula1
      HTML(paste(str1))
  })
  
  #HI formula 2
  output$global_het_indicies_formula_HI_2 = renderUI({
    myvals_obs = c()
    myvals_n = c()
    for(i in 1:input$num_pops){
      #observed heterozygotes
      myvals_obs = c(myvals_obs,  eval(parse(text = (paste0("input$slider2",i," /input$max",i)))))
      #N
      myvals_n = c(myvals_n, eval(parse(text = (paste0("input$max",i)))))
      
    }
    #reformat observed heterozygotes
    myvals_obs_reformatted = format(round(myvals_obs, 3), nsmall = 3)
    #add "x" and "+" to myvals_n and num_formula2.2, respectively
    x_var = c(rep(" x ", input$num_pops))
    plus_var = c(rep(" + ", input$num_pops))
    #combine numerator values and operators
    num = c(rbind(myvals_n, x_var, myvals_obs_reformatted, plus_var))
    den = c(rbind(myvals_n, plus_var))
    #remove the plus sign at end of num vector
    num = num[-length(num)]
    den = den[-length(den)]
    #Calculate HI
    HI = sum((myvals_obs * myvals_n)) / sum(myvals_n)
    #construct formula to print
    formula2 = c(paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',HI),"</font></b>"," = "),"<font color=\"#FF0000\">",num,"</font>","/","<font color=\"#FF0000\">",den,"</font>")
    #construct string outputs
    HTML(paste(formula2))
  })
  
  
  #HS formula 1
  output$global_het_indicies_formula_HS_1 = renderUI({
    num_formula1 = c()
    den_formula1 = c()
    for(i in 1:input$num_pops){
      num_formula1 = c(num_formula1, paste0("N","<sub>",i,"</sub>"," x <i>H</i>","<sub>","Exp",i,"</sub>"), " + ")
      den_formula1 = c(den_formula1, paste0("N","<sub>",i,"</sub>"), " + ")
    }
    #remove the plus sign at end of num_formula1 vector
    num_formula1 = num_formula1[-length(num_formula1)]
    den_formula1 = den_formula1[-length(den_formula1)]
    formula1 = c(paste0("<b><i>H</i><sub>S</sub></b>"," = "),num_formula1, "/", den_formula1)
    #construct string outputs
    str1 = formula1
    HTML(paste(str1))
  })
  
  #HS formula 2
  output$global_het_indicies_formula_HS_2 = renderUI({
    myvals_obs = c()
    myvals_exp = c()
    myvals_n = c()
    for(i in 1:input$num_pops){
      #expected heterozygotes
      myvals_exp = c(myvals_exp,  eval(parse(text = (paste0("2 * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")) * (1 - ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")))")))))
      #N
      myvals_n = c(myvals_n, eval(parse(text = (paste0("input$max",i)))))
    }
    #reformat myvals_exp
    myvals_exp_formatted =  format(round(myvals_exp, 3), nsmall = 3)
    #add "x" and "+" to myvals_n and myvals_exp, respectively
    x_var = c(rep(" x ", input$num_pops))
    plus_var = c(rep(" + ", input$num_pops))
    #combine numerator values and operators
    num = c(rbind(myvals_n, x_var, myvals_exp_formatted, plus_var))
    den = c(rbind(myvals_n, plus_var))
    #remove the plus sign at end of num vector
    num = num[-length(num)]
    den = den[-length(den)]
    #Calculate HS
    HS = sum((myvals_exp * myvals_n)) / sum(myvals_n)
    #construct formula to print
    formula2 = c(paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',HS),"</font></b>"," = "),"<font color=\"#FF0000\">",num,"</font>","/","<font color=\"#FF0000\">",den,"</font>")
    #construct string outputs
    HTML(paste(formula2))
  })
  
  #HT formula 1
  output$global_het_indicies_formula_HT_1 = renderUI({
    formula1 = c(paste0("<b><i>H</i><sub>T</sub></b>"," = 2 x pbar x qbar"))
    #construct string outputs
    HTML(paste(formula1))
  })
  #HT formula 1
  output$global_het_indicies_formula_HT_2 = renderUI({
    num = c()
    den = c()
    for(i in 1:input$num_pops){
      num = c(num, eval(parse(text = (paste0("input$max",i," * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,"))")))))
      den = c(den, eval(parse(text = (paste0("input$max",i)))))
    }
    #Calculate pbar and qbar
    pbar = sum(num) / sum(den)
    qbar = 1 - pbar
    #Calculate HT
    HT = 2*pbar*qbar
    #construct formula to print
    formula2 = c(paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',HT),"</b>"," = 2 x ",sprintf(fmt = '%#.3f',pbar)," x ",sprintf(fmt = '%#.3f',qbar),"</font>"))
    #construct string outputs
    HTML(paste(formula2))
  })
  
  #heterozgosity index table
  output$global_het_indicies_table = renderTable({
    m = matrix(nrow = 3, ncol = 1)
    myvals_obs = c()
    myvals_exp = c()
    myvals_n = c()
    num = c()
    den = c()
    for(i in 1:input$num_pops){
      myvals_obs = c(myvals_obs,  eval(parse(text = (paste0("input$slider2",i," /input$max",i)))))
      myvals_exp = c(myvals_exp,  eval(parse(text = (paste0("2 * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")) * (1 - ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")))")))))
      myvals_n = c(myvals_n, eval(parse(text = (paste0("input$max",i)))))
      num = c(num, eval(parse(text = (paste0("input$max",i," * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,"))")))))
      den = c(den, eval(parse(text = (paste0("input$max",i)))))
    }
    pbar = sum(num) / sum(den)
    #Hi
    m[1,1] = sum((myvals_obs * myvals_n)) / sum(myvals_n)
    #Hs
    m[2,1] = sum((myvals_exp * myvals_n)) / sum(myvals_n)
    #Ht
    m[3,1] = 2 * pbar * (1-pbar)
    colnames(m) = "values"
    m = as.data.frame(m)
    d = data.frame(heterozygosity_index = c("Hi", "Hs", "Ht"))
    df = cbind(d,m)
    df
  })#end of global het indicies table
  
  #Step 6.  CALCULATE THE GLOBAL  F-STATISTICS:
  # Compare and contrast the global Fis below with the local inbreeding coefficient Fs of Step 5.
  # Here we are using a weighted average of the individual heterozygosities over all the subpopulations.
  # Both FIS and  Fs are based on the observed heterozygosities, whereas Fst and Fit are based on expected heterozygosities.
  # Fis = (Hs - Hi) / Hs
  # Fst = (Ht - Hs) / Ht
  # Fit = (Ht - Hi) / Ht
  
  #Fstats formula
  output$global_fstats_formula = renderUI({
    myvals_obs = c()
    myvals_exp = c()
    myvals_n = c()
    num = c()
    den = c()
    for(i in 1:input$num_pops){
      myvals_obs = c(myvals_obs,  eval(parse(text = (paste0("input$slider2",i," /input$max",i)))))
      myvals_exp = c(myvals_exp,  eval(parse(text = (paste0("2 * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")) * (1 - ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")))")))))
      myvals_n = c(myvals_n, eval(parse(text = (paste0("input$max",i)))))
      num = c(num, eval(parse(text = (paste0("input$max",i," * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,"))")))))
      den = c(den, eval(parse(text = (paste0("input$max",i)))))
    }
    pbar = sum(num) / sum(den)
    #Hi
    HI = sum((myvals_obs * myvals_n)) / sum(myvals_n)
    #Hs
    HS = sum((myvals_exp * myvals_n)) / sum(myvals_n)
    #Ht
    HT = 2 * pbar * (1-pbar)
    #Construct formula to print
    #Fis_formula1
    str1 = paste0("<b><i>F</i><sub>IS</sub></b>"," = <i>H</i><sub>S</sub> - <i>H</i><sub>I</sub> / <i>H</i><sub>S</sub>")
    #Fis_formula2
    str2 = c(paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',(HS - HI)/HS),"</b>"," = ",sprintf(fmt = '%#.3f',HS)," - ",sprintf(fmt = '%#.3f',HI)," / ",sprintf(fmt = '%#.3f',HS) ,"</font>"))
    
    #Fst_formula1
    str3 = paste0("<b><i>F</i><sub>ST</sub></b>"," = <i>H</i><sub>T</sub> - <i>H</i><sub>S</sub> / <i>H</i><sub>T</sub>")
    #Fst_formula2
    str4 = c(paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',(HT - HS)/HT),"</b>"," = ",sprintf(fmt = '%#.3f',HT)," - ",sprintf(fmt = '%#.3f',HS)," / ",sprintf(fmt = '%#.3f',HT) ,"</font>"))
    
    #Fit_formula1
    str5 = paste0("<b><i>F</i><sub>IT</sub></b>"," = <i>H</i><sub>T</sub> - <i>H</i><sub>I</sub> / <i>H</i><sub>T</sub>")
    #Fit_formula2
    str6 = c(paste0("<font color=\"#FF0000\"><b>",sprintf(fmt = '%#.3f',(HT - HI)/HT),"</b>"," = ",sprintf(fmt = '%#.3f',HT)," - ",sprintf(fmt = '%#.3f',HI)," / ",sprintf(fmt = '%#.3f',HT) ,"</font>"))
  
    HTML(paste(str1, str2, str3, str4, str5, str6, sep = '</br>'))
  })
  
  output$global_fstats_table = renderTable({
    m = matrix(nrow = 3, ncol = 1)
    myvals_obs = c()
    myvals_exp = c()
    myvals_n = c()
    num = c()
    den = c()
    for(i in 1:input$num_pops){
      myvals_obs = c(myvals_obs,  eval(parse(text = (paste0("input$slider2",i," /input$max",i)))))
      myvals_exp = c(myvals_exp,  eval(parse(text = (paste0("2 * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")) * (1 - ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,")))")))))
      myvals_n = c(myvals_n, eval(parse(text = (paste0("input$max",i)))))
      num = c(num, eval(parse(text = (paste0("input$max",i," * ((2*input$slider1",i,"+input$slider2",i,")/(2*input$max",i,"))")))))
      den = c(den, eval(parse(text = (paste0("input$max",i)))))
    }
    pbar = sum(num) / sum(den)
    #Hi
    hi = sum((myvals_obs * myvals_n)) / sum(myvals_n)
    #Hs
    hs = sum((myvals_exp * myvals_n)) / sum(myvals_n)
    #Ht
    ht = 2 * pbar * (1-pbar)
    
    #Fis
    m[1,1] = (hs - hi) / hs
    #Fst
    m[2,1] = (ht - hs) / ht
    #Fit
    m[3,1] = (ht - hi) / ht
    colnames(m) = "values"
    m = as.data.frame(m)
    d = data.frame(fstat = c("Fis", "Fst", "Fit"))
    df = cbind(d,m)
    df
  })#end of global fstats table

  
  
}#end of server

shinyApp(ui = ui, server = server)
