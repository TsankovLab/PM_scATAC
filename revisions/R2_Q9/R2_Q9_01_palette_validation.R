# R2_Q9 -- old vs new sample palette: swatches under normal vision and simulated colour blindness, and the
# smallest pairwise colour distances (OKLab x 100; below ~15 normal / ~8 colour-blind two colours are easily confused)
suppressPackageStartupMessages({ library(farver); library(colorspace) })
source ('R2_Q9_palette_sample.R')
lab = function (hex) convert_colour (t(col2rgb (hex)), 'rgb', 'oklab')
pairD = function (hex)
  {
  ij = t(combn (length (hex), 2))
  f = function (h) { m = lab (h); sqrt (rowSums ((m[ij[,1],,drop=F] - m[ij[,2],,drop=F])^2)) * 100 }
  data.frame (a = names (hex)[ij[,1]], b = names (hex)[ij[,2]], normal = f (hex), deutan = f (deutan (hex)), protan = f (protan (hex)))
  }
tum = sample_scS_order
sets = list (old_tumours = palette_sample_old[tum], new_tumours = palette_sample[tum],
  old_all = palette_sample_old, new_all = palette_sample[names (palette_sample) != 'normal_pleura'])
res = do.call (rbind, lapply (names (sets), function(n) { d = pairD (sets[[n]]); cvd = pmin (d$deutan, d$protan)
  data.frame (set = n, n_colours = length (sets[[n]]), min_normal = round (min (d$normal), 1), worst_normal = paste (d$a, d$b)[which.min (d$normal)],
    pairs_normal_below_15 = sum (d$normal < 15), min_colourblind = round (min (cvd), 1), worst_colourblind = paste (d$a, d$b)[which.min (cvd)], pairs_colourblind_below_8 = sum (cvd < 8)) }))
print (res)
write.csv (res, 'R2_Q9_palette_distances.csv', row.names = FALSE)

pdf (file.path ('Plots', 'R2_Q9_palette_old_vs_new.pdf'), width = 8.5, height = 4.2)
ord = c(tum, 'P11_HOX', 'normal1', 'normal2', 'normal3', 'HU37', 'HU62')
par (mar = c(1, 9, 3, 1), xpd = NA)
rows = list ('old' = palette_sample_old[ord], 'old, deuteranopia' = deutan (palette_sample_old[ord]), 'old, protanopia' = protan (palette_sample_old[ord]),
  'new' = palette_sample[ord], 'new, deuteranopia' = deutan (palette_sample[ord]), 'new, protanopia' = protan (palette_sample[ord]))
plot (NA, xlim = c(0, length (ord)), ylim = c(0, 6.6), axes = FALSE, xlab = '', ylab = '')
y = c(5.6, 4.75, 3.9, 2.4, 1.55, 0.7)
for (i in seq_along (rows)) { rect (seq_along (ord) - 1 + .04, y[i], seq_along (ord) - .04, y[i] + .75, col = rows[[i]], border = NA)
  text (-.2, y[i] + .37, names (rows)[i], adj = 1, cex = .85, font = ifelse (i %in% c(1, 4), 2, 1)) }
text (seq_along (ord) - .5, 6.55, ord, cex = .7, srt = 45, adj = 0)
segments (c(0, 7, 12), 3.5, c(6.9, 10.9, 16.9), 3.5, lwd = 1.5, col = 'grey30')
text (c(3.45, 8.95, 14.45), 3.3, c('epithelioid-like', 'higher scS', 'subclone / non-tumour'), cex = .7, col = 'grey30')
dev.off()
