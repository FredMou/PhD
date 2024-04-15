

setwd('/Users/zhuofanmou/Documents/Phd/papers/paper2/barplot_events/TCGA')
#setwd('/Users/zhuofanmou/Documents/Phd/papers/paper2/barplot_events/ClariomD')
getwd()


stringsAsFactors = FALSE
stringsAsFactors = FALSE
# libraries:
library(readxl)
library(dplyr)

########################################################################################################
                                 # significant AS events from TCGA SpliceSeq #
########################################################################################################

# read in data
TCGA_spliceseq_DEAS <- read.delim("~/Documents/Phd/papers/paper2/barplot_events/TCGA/DEAS.txt")
df_sig <- as.data.frame(TCGA_spliceseq_DEAS)
class(df_sig)
View(df_sig)   

df_sig <- dplyr::filter(df_sig,  !is.na(gene)) # remove row with NA gene symbols

# # Summarisation of unique genes
# table(df_sig$`Event Type`)
# 
# # total
# gene_total <- unique(df_sig$`Gene name`)
# length(gene_total)
# 
# # ME
# gene_ME <- df_sig[ df_sig$`Event Type` =="Mutually Exclusive Exons", ]
# View(gene_ME)
# gene_ME <- unique(gene_ME$`Gene name`)
# length(gene_ME)
# 
# # A3SS
# gene_A3SS <- df_sig[ df_sig$`Event Type` =="Alternative 3' Splice Site", ]
# View(gene_A3SS)
# gene_A3SS <- unique(gene_A3SS$`Gene name`)
# length(gene_A3SS)
# 
# # A5SS
# gene_A5SS <- df_sig[ df_sig$`Event Type` =="Alternative 5' Splice Site", ]
# View(gene_A5SS)
# gene_A5SS <- unique(gene_A5SS$`Gene name`)
# length(gene_A5SS)
# 
# # CE
# gene_CE <- df_sig[ df_sig$`Event Type` =="Cassette Exon", ]
# View(gene_CE)
# gene_CE <- unique(gene_CE$`Gene name`)
# length(gene_CE)
# 
# # ALE
# gene_ALE <- df_sig[ df_sig$`Event Type` =="Alternative Last Exon", ]
# View(gene_ALE)
# gene_ALE <- unique(gene_ALE$`Gene name`)
# length(gene_ALE)
# 
# #AFE
# gene_AFE <- df_sig[ df_sig$`Event Type` =="Alternative First Exon", ]
# View(gene_AFE)
# gene_AFE <- unique(gene_AFE$`Gene name`)
# length(gene_AFE)
# 
# # RI
# gene_RI <- df_sig[ df_sig$`Event Type` =="Retained Intron", ]
# View(gene_RI)
# gene_RI <- unique(gene_RI$`Gene name`)
# length(gene_RI)
# 
# # complex event
# gene_ComplexEve <- df_sig[ df_sig$`Event Type` =="Complex Event", ]
# View(gene_ComplexEve)
# gene_ComplexEve <- unique(gene_ComplexEve$`Gene name`)
# length(gene_ComplexEve)
  
  
# table the event and corresponding number
table(df_sig$asType)
plot_df_sig <- as.data.frame(table(df_sig$asType))
View(plot_df_sig)


plot_df_sig$percent <- 100*(plot_df_sig$Freq/sum(plot_df_sig$Freq)) # Percentage of events of each type/all events

colnames(plot_df_sig) <- c('AS_Event', 'n', 'Percent') # rename the colnames


plt <- ggplot(plot_df_sig)+
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:5) * 500), # grid range
    color = "lightgrey"
  ) + 
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(str_wrap(AS_Event, 5), n),
      y = n,
      fill = Percent
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = .9
  ) +
  # # Add dots to represent the mean gain
  # geom_point(
  #   aes(
  #     x = reorder(str_wrap(region, 5),sum_length),
  #     y = mean_gain
  #   ),
  #   size = 3,
  #   color = "gray12"
  # ) +
  
  # Lollipop shaft for mean gain per region
  geom_segment(
    aes(
      x = reorder(str_wrap(AS_Event, 5), n),
      y = 0,
      xend = reorder(str_wrap(AS_Event, 5), n),
      yend = 2500
    ),
    linetype = "dashed",
    color = "gray12"
  ) +

  # Make it circular!
  coord_polar()

plt



plt <- plt +
  # Annotate the bars and the lollipops so the reader understands the scaling
  annotate(
    x = 9, 
    y = 100,
    label = "Number of AS Event",
    geom = "text",
    angle = -67.5,
    color = "gray12",
    size = 3.5#,
    # family = "Bell MT"
  ) +
  # annotate(
  #   x = 2,
  #   y = 3150,
  #   label = "Cummulative Length [FT]",
  #   geom = "text",
  #   angle = 23,
  #   color = "gray12",
  #   size = 2.5#,
  #   #family = "Bell MT"
  # ) +
  
  # # Annotate custom scale inside plot
  # annotate(
  #   x = 2, 
  #   y = 1100, 
  #   label = "1000", 
  #   geom = "text", 
  #   color = "gray12"#, 
  #   #family = "Bell MT"
  # ) +
  # annotate(
  #   x = 2, 
  #   y = 2100, 
  #   label = "2000", 
  #   geom = "text", 
  #   color = "gray12"#, 
  #   #family = "Bell MT"
  # ) +
  # annotate(
  #   x = 2, 
  #   y =3100, 
  #   label = "3000", 
  #   geom = "text", 
  #   color = "gray12"#, 
  #   #family = "Bell MT"
  # ) +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-800, 2500),
    expand = c(0, 0),
    breaks = c(0,500,1000,1500,2000,2500)
  ) +
  # New fill and legend title for number of tracks per region
  scale_fill_gradientn(
    "% of AS Events",
    colours = c( "#492a21","#a14e48","#e61942","#ffae34")#c( "#5E3023","#C08552","#DAB49D","#F3E9DC") # c( "#492a21","#a14e48","#e61942","#ffae34") # c( "#6C5B7B","#C06C84","#F67280","#F8B195")
  ) +
  # # Make the guide for the fill discrete
  # guides(
  #   fill = guide_colorsteps(
  #     barwidth = 15, barheight = .5, title.position = "top", title.hjust = .5
  #   )
  # ) +
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 12),
    # Move the legend to the bottom
    legend.position = "bottom",
  )

plt






plt <- plt +
  # Add labels
  labs(
    title = "\nSignificant Events-TCGA-PRAD Dataset") +
    # subtitle = paste(
    #   "\nThis Visualisation shows the cummulative length of tracks,",
    #   "the amount of tracks and the mean gain in elevation per location.\n",
    #   "If you are an experienced hiker, you might want to go",
    #   "to the North Cascades since there are a lot of tracks,",
    #   "higher elevations and total length to overcome.",
    #   sep = "\n"
    # ),
    # caption = "\n\nData Visualisation by Tobias Stalder\ntobias-stalder.netlify.app\nSource: TidyX Crew (Ellis Hughes, Patrick Ward)\nLink to Data: github.com/rfordatascience/tidytuesday/blob/master/data/2020/2020-11-24/readme.md") +
  # Customize general theme
  theme(

    # Set default color and font family for the text
    text = element_text(color = "gray12"#,
                        # family = "Bell MT"
    ),

    # Customize the text in the title, subtitle, and caption
    plot.title = element_text(face = "bold", size = 25, hjust = 0.05),
    plot.subtitle = element_text(size = 14, hjust = 0.05),
    plot.caption = element_text(size = 10, hjust = .5),

    # Make the background white and remove extra grid lines
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_blank()
  )
# Use `ggsave("plot.png", plt,width=9, height=12.6)` to save it as in the output
plt

ggsave("TCGAspliceseq_8440_sigEvents_labelledGen_contLegend.pdf", plt,width=9, height=12.6)






########################################################################################################
                              # DS parent genes of significant events #
########################################################################################################

# read in data
TCGA_spliceseq_DEAS <- read.delim("~/Documents/Phd/papers/paper2/barplot_events/TCGA/DEAS.txt")
df_sig <- as.data.frame(TCGA_spliceseq_DEAS)
class(df_sig)
View(df_sig)   

df_sig <- dplyr::filter(df_sig,  !is.na(gene)) # remove row with NA gene symbols; 1849
length(unique(df_sig$gene)) # 4257


df_sig <- df_sig %>%
       distinct(asType, gene) %>%
       group_by(asType) %>%
       summarise("Event in Genes" = n())
View(df_sig)


# table the event and corresponding number
#table(df_sig$`Event Type`)
plot_df_sig <- df_sig #as.data.frame(table(df_sig$`Event Type`))
View(plot_df_sig)


plot_df_sig$percent <- 100*(plot_df_sig$`Event in Genes`/sum(plot_df_sig$`Event in Genes`)) # Percentage of events of each type/all events

colnames(plot_df_sig) <- c('AS_Event', 'n', 'Percent') # rename the colnames


plt <- ggplot(plot_df_sig)+
  # Make custom panel grid
  geom_hline(
    aes(yintercept = y), 
    data.frame(y = c(0:8) * 200), # grid range
    color = "lightgrey"
  ) + 
  # Add bars to represent the cumulative track lengths
  # str_wrap(region, 5) wraps the text so each line has at most 5 characters
  # (but it doesn't break long words!)
  geom_col(
    aes(
      x = reorder(str_wrap(AS_Event, 5), n),
      y = n,
      fill = Percent
    ),
    position = "dodge2",
    show.legend = TRUE,
    alpha = .9
  ) +
  # # Add dots to represent the mean gain
  # geom_point(
  #   aes(
  #     x = reorder(str_wrap(region, 5),sum_length),
  #     y = mean_gain
  #   ),
  #   size = 3,
  #   color = "gray12"
  # ) +
  
  # Lollipop shaft for mean gain per region
geom_segment(
  aes(
    x = reorder(str_wrap(AS_Event, 5), n),
    y = 0,
    xend = reorder(str_wrap(AS_Event, 5), n),
    yend = 1600
  ),
  linetype = "dashed",
  color = "gray12"
) +
  
  # Make it circular!
  coord_polar()

plt



plt <- plt +
  # Annotate the bars and the lollipops so the reader understands the scaling
  annotate(
    x = 9, 
    y = 100,
    label = "Number of AS Event",
    geom = "text",
    angle = -67.5,
    color = "gray12",
    size = 3.5#,
    # family = "Bell MT"
  ) +
  # annotate(
  #   x = 2,
  #   y = 3150,
  #   label = "Cummulative Length [FT]",
  #   geom = "text",
  #   angle = 23,
  #   color = "gray12",
  #   size = 2.5#,
  #   #family = "Bell MT"
  # ) +
  
# # Annotate custom scale inside plot
# annotate(
#   x = 2, 
#   y = 1100, 
#   label = "1000", 
#   geom = "text", 
#   color = "gray12"#, 
#   #family = "Bell MT"
# ) +
# annotate(
#   x = 2, 
#   y = 2100, 
#   label = "2000", 
#   geom = "text", 
#   color = "gray12"#, 
#   #family = "Bell MT"
# ) +
# annotate(
#   x = 2, 
#   y =3100, 
#   label = "3000", 
#   geom = "text", 
#   color = "gray12"#, 
#   #family = "Bell MT"
# ) +
# Scale y axis so bars don't start in the center
scale_y_continuous(
  limits = c(-500, 1600),
  expand = c(0, 0),
  breaks = c(0,200,400,600,800,1000,1200,1400,1600)
) +
  # New fill and legend title for number of tracks per region
  scale_fill_gradientn(
    "% of AS Events",
    colours = c( "#492a21","#a14e48","#e61942","#ffae34")#c( "#5E3023","#C08552","#DAB49D","#F3E9DC") # c( "#492a21","#a14e48","#e61942","#ffae34") # c( "#6C5B7B","#C06C84","#F67280","#F8B195")
  ) +
  # # Make the guide for the fill discrete
  # guides(
  #   fill = guide_colorsteps(
  #     barwidth = 15, barheight = .5, title.position = "top", title.hjust = .5
  #   )
  # ) +
  theme(
    # Remove axis ticks and text
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    # Use gray text for the region names
    axis.text.x = element_text(color = "gray12", size = 12),
    # Move the legend to the bottom
    legend.position = "bottom",
  )

plt






plt <- plt +
  # Add labels
  labs(
    title = "\nSignificant DSG by Event Type-Clariom D Dataset") +
  # subtitle = paste(
  #   "\nThis Visualisation shows the cummulative length of tracks,",
  #   "the amount of tracks and the mean gain in elevation per location.\n",
  #   "If you are an experienced hiker, you might want to go",
  #   "to the North Cascades since there are a lot of tracks,",
  #   "higher elevations and total length to overcome.",
  #   sep = "\n"
  # ),
  # caption = "\n\nData Visualisation by Tobias Stalder\ntobias-stalder.netlify.app\nSource: TidyX Crew (Ellis Hughes, Patrick Ward)\nLink to Data: github.com/rfordatascience/tidytuesday/blob/master/data/2020/2020-11-24/readme.md") +
  # Customize general theme
  theme(
    
    # Set default color and font family for the text
    text = element_text(color = "gray12"#,
                        # family = "Bell MT"
    ),
    
    # Customize the text in the title, subtitle, and caption
    plot.title = element_text(face = "bold", size = 25, hjust = 0.05),
    plot.subtitle = element_text(size = 14, hjust = 0.05),
    plot.caption = element_text(size = 10, hjust = .5),
    
    # Make the background white and remove extra grid lines
    panel.background = element_rect(fill = "white", color = "white"),
    panel.grid = element_blank(),
    panel.grid.major.x = element_blank()
  )
# Use `ggsave("plot.png", plt,width=9, height=12.6)` to save it as in the output
plt

ggsave("TCGAspliceseq_5428_sigDSG_labelledGen_contLegend.pdf", plt,width=9, height=12.6)

