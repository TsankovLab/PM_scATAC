library (paletteer)
library (circlize)

# Load palettes 
palette_chromvar = paletteer::paletteer_c("pals::kovesi.diverging_bwr_40_95_c42",100)
palette_chromvar_fun = colorRamp2(c(-20,0,20), c(palette_chromvar[1],'white',palette_chromvar[100]))
