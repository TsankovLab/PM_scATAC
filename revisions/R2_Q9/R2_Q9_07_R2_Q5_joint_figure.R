# R2_Q9 -- joint scATAC/scRNA embedding figure of R2_Q5 (09_figure_joint.R) with the improved sample palette.
# The R2_Q5 script is run unchanged except for the palette and the output folder.
q9  = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q9'
q5  = '/sc/arion/projects/Tsankov_Normal_Lung/Bruno/mesothelioma/scATAC_PM/git_repo_claude/R2_Q5'
setwd (q5)
code = readLines ('09_figure_joint.R')
i = grep ('^source\\("00_common.R"\\)', code); stopifnot (length (i) == 1)
code = append (code, sprintf ('source ("%s")', file.path (q9, 'R2_Q9_palette_sample.R')), after = i)
stopifnot (any (grepl ('"Plots/R2Q5_joint', code, fixed = TRUE)))
code = gsub ('"Plots/R2Q5_joint', sprintf ('"%s/R2Q5_joint', file.path (q9, 'Plots')), code, fixed = TRUE)
eval (parse (text = code))
