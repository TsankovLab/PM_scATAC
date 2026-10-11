# R2_Q9 -- improved sample palette (drop-in replacement of palette_sample in git_repo/utils/palettes.R)
#
# Logic, so that the colour tells something at a glance:
#   * tumours: 11 clearly different hues (chosen in R2_Q9_00_palette_design.R; smallest pairwise OKLab
#     distance 18 for normal vision and 10 under simulated red-green colour blindness, against 5 and 2 for
#     the old palette). They are assigned along the sarcomatoid score (median scS of the malignant cells,
#     scATAC): cool colours (green, teal, blue) = epithelioid-like tumours, violet/purple = intermediate,
#     warm colours (pink, amber, red) = tumours with higher scS.
#   * non-tumour samples: neutral black/greys, so that "normal" is recognised without reading the legend, each
#     with its own shade (the old palette gave the same red to all of them). Only near-black stays
#     clear of every tumour colour under colour-blind simulation, so normal1 (the normal sample shown together
#     with the tumours) is near-black; the mid greys of the other normal samples can come close to P10/P12.
#   * P11_HOX (HOX-high subclone of P11): a light tint of the P11 colour.
palette_sample = c(
  P8  = '#1C9B00',  # green        most epithelioid-like
  P4  = '#00E58D',  # mint
  P10 = '#006557',  # dark teal
  P12 = '#469CB5',  # steel blue
  P23 = '#0057D7',  # blue
  P11 = '#926FFF',  # violet
  P5  = '#7E2491',  # purple
  P3  = '#FF84FD',  # pink
  P14 = '#FFA900',  # amber
  P13 = '#FF5866',  # coral
  P1  = '#AC0000',  # dark red     highest scS
  normal1 = '#222222', normal_pleura = '#222222',   # near-black: the normal sample plotted with the tumours
  normal2 = '#737373', normal3 = '#9A9A9A',
  HU37 = '#4A4A4A', HU62 = '#B8B8B8',
  P11_HOX = '#D8CCFF')
# same order of names as the old palette first, so that positional uses keep pointing at the same samples
palette_sample = palette_sample[c('P1','P13','P3','P12','P5','P11','P4','P8','P10','P14','P23','HU62','HU37','normal1','normal2','normal3','P11_HOX','normal_pleura')]

palette_sample_old = setNames (
  c('#4B2D4B','#3C3C5A','#4B6987','#789696','#968787','#D2C3C3','#D2B487','#875A2D','#5A3C2D','#873C3C','#00008B','#FF0000','#FF0000','#FF0000','#FF0000','#FF0000','#EE82EE'),
  c('P1','P13','P3','P12','P5','P11','P4','P8','P10','P14','P23','HU62','HU37','normal1','normal2','normal3','P11_HOX'))
# tumours ordered by median scS (scATAC, R2_Q10)
sample_scS_order = c('P8','P4','P10','P12','P23','P11','P5','P3','P14','P13','P1')
