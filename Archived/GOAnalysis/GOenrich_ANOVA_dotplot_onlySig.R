#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
#BiocManager::install('______')

#Load libraries
library(ggplot2)
library(gridExtra)

#Read in data on GO terms (BP, MF, and CC) for each ANOVA category 
UVvsVis_BP_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_UVvsVIS_output/2WayANOVA_UVvsVIS_all_genes_BP_GO_terms.csv', row.names = 1)
UVvsVis_MF_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_UVvsVIS_output/2WayANOVA_UVvsVIS_all_genes_MF_GO_terms.csv', row.names = 1)
UVvsVis_CC_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_UVvsVIS_output/2WayANOVA_UVvsVIS_all_genes_CC_GO_terms.csv', row.names = 1)

TvsN_BP_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_TvsN_output/2WayANOVA_TvsN_all_genes_BP_GO_terms.csv', row.names = 1)
TvsN_MF_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_TvsN_output/2WayANOVA_TvsN_all_genes_MF_GO_terms.csv', row.names = 1)
TvsN_CC_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_TvsN_output/2WayANOVA_TvsN_all_genes_CC_GO_terms.csv', row.names = 1)

interaction_BP_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_interaction_output/2WayANOVA_interaction_all_genes_BP_GO_terms.csv', row.names = 1)
interaction_MF_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_interaction_output/2WayANOVA_interaction_all_genes_MF_GO_terms.csv', row.names = 1)
interaction_CC_GO_terms <- read.csv('/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_results/GOenrich_interaction_output/2WayANOVA_interaction_all_genes_CC_GO_terms.csv', row.names = 1)

#--------------------------------------------------------------------------------------------------

#ALL BP WORK
    #Filter for top 5 significant in each (get the 1st 5 rows...then take only the significant ones).
    #They should already be in decreasing significance order.
UVvsVis_BP_GO_top5 <- UVvsVis_BP_GO_terms[1:5, ]
TvsN_BP_GO_top5 <- TvsN_BP_GO_terms[1:5, ]
interaction_BP_GO_top5 <- interaction_BP_GO_terms[1:5, ]

UVvsVis_BP_GO_sig <- UVvsVis_BP_GO_top5[which(UVvsVis_BP_GO_top5$p_adjusted <= 0.05), ]
TvsN_BP_GO_sig <- TvsN_BP_GO_top5[which(TvsN_BP_GO_top5$p_adjusted <= 0.05), ]
interaction_BP_GO_sig <- interaction_BP_GO_top5[which(interaction_BP_GO_top5$p_adjusted <= 0.05), ]


    #List of all 30 BP for each effect (duplicates removed)
list_all_BP_GO_included <- unique(c(UVvsVis_BP_GO_sig$GO.ID, TvsN_BP_GO_sig$GO.ID, interaction_BP_GO_sig$GO.ID))


    #Add a column labeling the effect to each 
#UVvsVis_BP_plot_table <- cbind("Effect" = 'Treatment', UVvsVis_BP_GO_sig)
#TvsN_BP_plot_table <- cbind("Effect" = 'Tolerance', TvsN_BP_GO_sig)
interaction_BP_plot_table <- cbind("Effect" = 'Interaction', interaction_BP_GO_sig)

#--------------------------------------------------------------------------------------------------

#ALL MF WORK
    #Filter for top 5 significant in each (get the 1st 5 rows...then take only the significant ones).
    #They should already be in decreasing significance order.
UVvsVis_MF_GO_top5 <- UVvsVis_MF_GO_terms[1:5, ]
TvsN_MF_GO_top5 <- TvsN_MF_GO_terms[1:5, ]
interaction_MF_GO_top5 <- interaction_MF_GO_terms[1:5, ]

UVvsVis_MF_GO_sig <- UVvsVis_MF_GO_top5[which(UVvsVis_MF_GO_top5$p_adjusted <= 0.05), ]
TvsN_MF_GO_sig <- TvsN_MF_GO_top5[which(TvsN_MF_GO_top5$p_adjusted <= 0.05), ]
interaction_MF_GO_sig <- interaction_MF_GO_top5[which(interaction_MF_GO_top5$p_adjusted <= 0.05), ]


    #List of all 30 MF for each effect (duplicates removed)
list_all_MF_GO_included <- unique(c(UVvsVis_MF_GO_sig$GO.ID, TvsN_MF_GO_sig$GO.ID, interaction_MF_GO_sig$GO.ID))


    #Add a column labeling the effect to each 
UVvsVis_MF_plot_table <- cbind("Effect" = 'Treatment', UVvsVis_MF_GO_sig)
TvsN_MF_plot_table <- cbind("Effect" = 'Tolerance', TvsN_MF_GO_sig)
interaction_MF_plot_table <- cbind("Effect" = 'Interaction', interaction_MF_GO_sig)

#--------------------------------------------------------------------------------------------------

#ALL CC WORK
    #Filter for top 5 significant in each (get the 1st 5 rows...then take only the significant ones).
    #They should already be in decreasing significance order.
UVvsVis_CC_GO_top5 <- UVvsVis_CC_GO_terms[1:5, ]
TvsN_CC_GO_top5 <- TvsN_CC_GO_terms[1:5, ]
interaction_CC_GO_top5 <- interaction_CC_GO_terms[1:5, ]


UVvsVis_CC_GO_sig <- UVvsVis_CC_GO_top5[which(UVvsVis_CC_GO_top5$p_adjusted <= 0.05), ]
TvsN_CC_GO_sig <- TvsN_CC_GO_top5[which(TvsN_CC_GO_top5$p_adjusted <= 0.05), ]
interaction_CC_GO_sig <- interaction_CC_GO_top5[which(interaction_CC_GO_top5$p_adjusted <= 0.05), ]

    #List of all 30 CC for each effect (duplicates removed)
list_all_CC_GO_included <- unique(c(UVvsVis_CC_GO_sig$GO.ID, TvsN_CC_GO_sig$GO.ID, interaction_CC_GO_sig$GO.ID))


    #Add a column labeling the effect to each 
UVvsVis_CC_plot_table <- cbind("Effect" = 'Treatment', UVvsVis_CC_GO_sig)
TvsN_CC_plot_table <- cbind("Effect" = 'Tolerance', TvsN_CC_GO_sig)
interaction_CC_plot_table <- cbind("Effect" = 'Interaction', interaction_CC_GO_sig)

#--------------------------------------------------------------------------------------------------

#Combine all tables into 1
all_effect_BP_plot_table <- interaction_BP_plot_table   # UVvsVis_MF_plot_table, TvsN_BP_plot_table weren't included because no significant terms were in the table 
all_effect_BP_plot_table <- cbind('GO_cat' = 'BP', all_effect_BP_plot_table)

all_effect_MF_plot_table <- rbind(UVvsVis_MF_plot_table, TvsN_MF_plot_table, interaction_MF_plot_table)
all_effect_MF_plot_table <- cbind('GO_cat' = 'MF', all_effect_MF_plot_table)

all_effect_CC_plot_table <- rbind(UVvsVis_CC_plot_table, TvsN_CC_plot_table, interaction_CC_plot_table)
all_effect_CC_plot_table <- cbind('GO_cat' = 'CC', all_effect_CC_plot_table)

all_effect_plot_table <- rbind(all_effect_BP_plot_table, all_effect_MF_plot_table, all_effect_CC_plot_table)

#--------------------------------------------------------------------------------------------------

#Make plot
x_axis_order <- factor(all_effect_plot_table$Effect, levels = c('Treatment', 'Tolerance', 'Interaction'))
facet <- factor(all_effect_plot_table$GO_cat, levels = c('BP', 'MF', 'CC'))

p <- ggplot(data = all_effect_plot_table, aes(x = x_axis_order, y = Term, size = Significant, color = p_adjusted)) 
p <- p + facet_grid(rows = facet, space = 'free_y', scales = 'free')
p <- p + geom_point() + scale_color_gradientn(colors = heat.colors(10), limits=c(0, 0.05)) + theme_bw()
p <- p + xlab('ANOVA Effect') + ylab('GO Term') 
p <- p + labs(color = 'Adjusted p-value', size = 'Gene rank')
#--------------------------------------------------------------------------------------------------

#Save to file
final_plot <- p
ggsave('GOenrich_ANOVA_Dotplot_OnlySig_CustomAnnotation_AllGenes.pdf', plot = final_plot, device = 'pdf')
file.move(paste(getwd(), '/', "GOenrich_ANOVA_Dotplot_OnlySig_CustomAnnotation_AllGenes.pdf", sep = ''), "/Users/bryanmichalek/Documents/Notre_Dame/Spring 2021/Pfrender/ANOVA_Results/", overwrite = TRUE)

