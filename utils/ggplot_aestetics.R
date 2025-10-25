# Set ggplot theme aestetics
gtheme = theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )
gtheme_no_rot = theme(
      #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )


gtheme_italic = theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face='italic'),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

gtheme_no_text = theme(
      axis.text.x = element_blank(),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_blank(),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

gtheme_text = theme(
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_blank(),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )

gtheme_italic_y = theme(
      axis.text.y = element_text(face='italic'),
      axis.line =element_line(colour = 'black', size = .1),
        axis.ticks = element_line(colour = "black", size = .1),
      panel.background = element_blank()#,
    #panel.border = element_blank(),
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank()
  )


# Boxplots / violins ####
bxp = geom_boxplot (
    linewidth = .2,
    width=.8,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) 
bxpv = geom_boxplot (
    linewidth = .2,
    width=.4,
    outlier.alpha = 0.2,
    outlier.size = 1,
     size=0.6, alpha=0.7
     ) 

vlp = geom_violin (trim=TRUE,size=2,
    width=1,
    scale='width',
    linewidth = .2, alpha=0.7)

