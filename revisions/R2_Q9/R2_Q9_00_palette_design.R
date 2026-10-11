# R2_Q9 -- design of an improved sample palette
#
# The current palette_sample (git_repo/utils/palettes.R) takes ten muted brown/grey/blue tones from an
# impressionist palette for the tumours and one identical red for every normal sample, so samples are hard
# to tell apart. Here 11 tumour colours are chosen by maximising the smallest pairwise colour difference,
# measured in OKLab under normal vision and under simulated deuteranopia and protanopia (Machado 2009,
# colorspace; see NORMAL_WEIGHT below), within a lightness/chroma band that reads well on a white page.
# The chosen colours are written to the console; the final hex values are fixed in R2_Q9_palette_sample.R.

suppressPackageStartupMessages({ library(farver); library(colorspace) })
set.seed(1)

dE = function (a, b) sqrt (rowSums ((a - b)^2)) * 100              # OKLab distance x 100
lab = function (hex) convert_colour (t(col2rgb (hex)), 'rgb', 'oklab')
views = function (hex) list (normal = lab (hex), deutan = lab (deutan (hex)), protan = lab (protan (hex)))
min_pair = function (hex)
  {
  v = views (hex); ij = t(combn (length (hex), 2))
  d = sapply (v, function(m) dE (m[ij[,1],,drop=F], m[ij[,2],,drop=F]))
  d = matrix (d, ncol = 3, dimnames = list (NULL, names (v)))
  w = which.min (apply (d, 1, min))
  list (min = min (d), normal = min (d[,'normal']), cvd = min (d[,c('deutan','protan')]), worst = hex[ij[w,]])
  }

# candidates: OKLCH grid, lightness 0.45-0.80, chroma 0.09-0.20, in sRGB gamut
g = expand.grid (l = seq (0.45, 0.80, 0.05), c = c(0.09, 0.12, 0.15, 0.18, 0.21), h = seq (0, 350, 10))
rgbm = convert_colour (as.matrix (g), 'oklch', 'rgb')
ok = rowSums (rgbm < 0 | rgbm > 255) == 0
cand = unique (encode_colour (rgbm[ok,], from = 'rgb'))
V = views (cand)
# Score of a pair: the smaller of (normal-vision distance / NORMAL_WEIGHT) and the colour-blind distances.
# Optimising the colour-blind distance alone returns only blues and oranges, which are safe but look like two
# groups to most readers; asking normal vision for NORMAL_WEIGHT times the colour-blind distance keeps the full
# range of hues (red, green, purple) and still spreads lightness so that pairs stay apart when hue is lost.
NORMAL_WEIGHT = as.numeric (Sys.getenv ('NORMAL_WEIGHT', '1.8'))
Dn = as.matrix (dist (V$normal)) * 100
Dc = pmin (as.matrix (dist (V$deutan)), as.matrix (dist (V$protan))) * 100
D = pmin (Dn / NORMAL_WEIGHT, Dc)
message (length (cand), ' candidate colours')

# farthest-point selection from many random starts, followed by pairwise swaps
n = 11
best = NULL; best_s = -1
for (it in 1:400)
  {
  s = sample (length (cand), 1)
  while (length (s) < n) s = c(s, which.max (apply (D[, s, drop = FALSE], 1, min) * !(seq_along (cand) %in% s)))
  repeat
    {
    improved = FALSE
    for (k in seq_len (n))
      {
      cur = min (D[s, s][upper.tri (diag (n))])
      alt = apply (D[, s[-k], drop = FALSE], 1, min); alt[s] = -1
      if (max (alt) > min (D[s[k], s[-k]]) && min (c(max (alt), D[s[-k], s[-k]][upper.tri (diag (n - 1))])) > cur) { s[k] = which.max (alt); improved = TRUE }
      }
    if (!improved) break
    }
  sc = min (D[s, s][upper.tri (diag (n))])
  if (sc > best_s) { best_s = sc; best = s }
  }
new = cand[best]
h = convert_colour (t(col2rgb (new)), 'rgb', 'oklch')
new = new[order (h[,'h'])]
cat ('new tumour colours (by hue):', new, '\n')
print (round (convert_colour (t(col2rgb (new)), 'rgb', 'oklch'), 2))
r = min_pair (new); cat (sprintf ('NEW  min pairwise dE: normal %.1f | colour-blind %.1f | worst pair %s %s\n', r$normal, r$cvd, r$worst[1], r$worst[2]))

old = c('#4B2D4B','#3C3C5A','#4B6987','#789696','#968787','#D2C3C3','#D2B487','#875A2D','#5A3C2D','#873C3C','#00008B')
r = min_pair (old); cat (sprintf ('OLD  min pairwise dE: normal %.1f | colour-blind %.1f | worst pair %s %s\n', r$normal, r$cvd, r$worst[1], r$worst[2]))
