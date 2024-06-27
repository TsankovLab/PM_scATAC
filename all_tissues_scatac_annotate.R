archp$celltype = 0

# Map Tsankov cell type annotation ####
prj = 'Tsankov_lung'
archp_barcodes = rownames (archp@cellColData)[archp$project == prj]
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[2])
archp_barcodes = data.frame (barcode = archp_barcodes, sample = archp$Sample[archp$project == prj])
archp_barcodes = split (archp_barcodes, archp_barcodes$sample)

ext_meta = meta[[3]]
ext_meta$barcode = sapply (rownames(ext_meta), function(x) unlist (strsplit (x, '\\#'))[2])
ext_meta$sample = sapply (rownames(ext_meta), function(x) unlist (strsplit (x, '\\#'))[1])
ext_meta = split (ext_meta, ext_meta$sample)

match_sample = sapply (archp_barcodes, function(x) sapply (ext_meta, function(y) sum (x$barcode %in% y$barcode)))
names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = lapply (names(archp_barcodes), function(x) 
  {
  archp_barcodes[[x]]$celltype = ext_meta[[x]]$cell_type[match(archp_barcodes[[x]]$barcode, ext_meta[[x]]$barcode)]
  archp_barcodes[[x]]
  })
#names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = do.call (rbind, archp_barcodes)
archp$celltype[match(rownames(archp_barcodes), rownames(archp@cellColData))] = archp_barcodes$celltype

# Map JShendure cell type annotation ####
prj = 'JShendure'
archp_barcodes = rownames (archp@cellColData)[archp$project == prj]
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[2])
archp_barcodes = data.frame (barcode = archp_barcodes, sample = archp$Sample[archp$project == prj])
archp_barcodes = split (archp_barcodes, archp_barcodes$sample)

ext_meta = meta[[2]]
ext_meta$barcode = rownames(ext_meta)
ext_meta$sample = ext_meta$sample_name
ext_meta = split (ext_meta, ext_meta$sample)

match_sample = sapply (archp_barcodes, function(x) sapply (ext_meta, function(y) sum (x$barcode %in% y$barcode)))
names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = archp_barcodes[!names(archp_barcodes) == '']
archp_barcodes2 = lapply (names(archp_barcodes), function(x) 
  {
  archp_barcodes[[x]]$celltype = ext_meta[[x]]$cell_type[match(archp_barcodes[[x]]$barcode, ext_meta[[x]]$barcode)]
  archp_barcodes[[x]]
  })
archp_barcodes2 = do.call (rbind, archp_barcodes2)
archp$celltype[match(rownames(archp_barcodes2), rownames(archp@cellColData))] = archp_barcodes2$celltype

# Map BingRen cell type annotation ####
prj = 'bingren_pan'
archp_barcodes = rownames (archp@cellColData)[archp$project == prj]
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[2])
archp_barcodes = data.frame (barcode = archp_barcodes, sample = archp$Sample[archp$project == prj])
archp_barcodes = split (archp_barcodes, archp_barcodes$sample)

ext_meta = meta[[1]]
ext_meta$barcode = sapply (ext_meta$cellID, function(x) unlist (strsplit (x, '\\+'))[2])
ext_meta$sample = ext_meta$sample
ext_meta = split (ext_meta, ext_meta$sample)

match_sample = sapply (archp_barcodes, function(x) sapply (ext_meta, function(y) sum (x$barcode %in% y$barcode)))
names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = archp_barcodes[!names(archp_barcodes) == '']
archp_barcodes2 = lapply (names(archp_barcodes), function(x) 
  {
  archp_barcodes[[x]]$celltype = ext_meta[[x]]$cell.type[match(archp_barcodes[[x]]$barcode, ext_meta[[x]]$barcode)]
  archp_barcodes[[x]]
  })
archp_barcodes2 = do.call (rbind, archp_barcodes2)
archp$celltype[match(rownames(archp_barcodes2), rownames(archp@cellColData))] = archp_barcodes2$celltype


# Map greenleaf Brain cell type annotation ####
prj = 'greenleaf_brain'
archp_barcodes = rownames (archp@cellColData)[archp$project == prj]
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[2])
archp_barcodes = data.frame (barcode = archp_barcodes, sample = archp$Sample[archp$project == prj])
archp_barcodes = split (archp_barcodes, archp_barcodes$sample)

celltype_mapping = read.csv ('metadata_external_data/greenleaf_brain_celltype_mapping.csv')
ext_meta = meta[[4]]
ext_meta$barcode = ext_meta$Cell.Barcode
ext_meta$sample = ext_meta$Sample.ID
celltype_gf = setNames (celltype_mapping$class, celltype_mapping$X...cluster)
ext_meta$celltype = celltype_gf[ext_meta$Iterative.LSI.Clusters]
ext_meta = split (ext_meta, ext_meta$sample)

match_sample = sapply (archp_barcodes, function(x) sapply (ext_meta, function(y) sum (x$barcode %in% y$barcode)))
names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = archp_barcodes[!names(archp_barcodes) == '']
archp_barcodes2 = lapply (names(archp_barcodes), function(x) 
  {
  archp_barcodes[[x]]$celltype = ext_meta[[x]]$celltype[match(archp_barcodes[[x]]$barcode, ext_meta[[x]]$barcode)]
  archp_barcodes[[x]]
  })
archp_barcodes2 = do.call (rbind, archp_barcodes2)
archp_barcodes2 = archp_barcodes2[archp_barcodes2$barcode %in% rownames (archp@cellColData),]
archp$celltype[match(rownames(archp_barcodes2), rownames(archp@cellColData))] = archp_barcodes2$celltype


# Map kidney cell type annotation ####
prj = 'yang_kidney'
archp_barcodes = rownames (archp@cellColData)[archp$project == prj]
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[2])
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '-'))[1])
archp_barcodes = data.frame (barcode = archp_barcodes, sample = archp$Sample[archp$project == prj])
archp_barcodes = split (archp_barcodes, archp_barcodes$sample)

ext_meta = meta[[5]]
ext_meta$barcode = sapply (rownames(ext_meta), function(x) unlist (strsplit (x, '-'))[1])
ext_meta$sample = ext_meta$batch
ext_meta = split (ext_meta, ext_meta$sample)

match_sample = sapply (archp_barcodes, function(x) sapply (ext_meta, function(y) sum (x$barcode %in% y$barcode)))
names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = archp_barcodes[!names(archp_barcodes) == '']
archp_barcodes2 = lapply (names(archp_barcodes), function(x) 
  {
  archp_barcodes[[x]]$celltype = ext_meta[[x]]$celltype[match(archp_barcodes[[x]]$barcode, ext_meta[[x]]$barcode)]
  archp_barcodes[[x]]$barcode = rownames(archp_barcodes[[x]])
  archp_barcodes[[x]]
  })
archp_barcodes2 = do.call (rbind, archp_barcodes2)
archp$celltype[match(archp_barcodes2$barcode, rownames(archp@cellColData))] = archp_barcodes2$celltype


# Map fetal lung cell type annotation ####
prj = 'rawlins_fetal_lung'
archp_barcodes = rownames (archp@cellColData)[archp$project == prj]
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[2])
archp_barcodes = gsub ('-1','',archp_barcodes)
#archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[1])
archp_barcodes = data.frame (barcode = archp_barcodes, sample = archp$Sample[archp$project == prj])
archp_barcodes = split (archp_barcodes, archp_barcodes$sample)

ext_meta = meta[[7]]
ext_meta$barcode = sapply (rownames(ext_meta), function(x) unlist (strsplit (x, '\\#'))[2])
ext_meta$sample = ext_meta$sample_name
ext_meta = split (ext_meta, ext_meta$sample)

match_sample = sapply (archp_barcodes, function(x) sapply (ext_meta, function(y) sum (x$barcode %in% y$barcode)))
names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = archp_barcodes[!names(archp_barcodes) == '']
archp_barcodes2 = lapply (names(archp_barcodes), function(x) 
  {
  archp_barcodes[[x]]$celltype = ext_meta[[x]]$cell_type[match(archp_barcodes[[x]]$barcode, ext_meta[[x]]$barcode)]
  archp_barcodes[[x]]$barcode = rownames(archp_barcodes[[x]])
  archp_barcodes[[x]]
  })
archp_barcodes2 = do.call (rbind, archp_barcodes2)
archp$celltype[match(archp_barcodes2$barcode, rownames(archp@cellColData))] = archp_barcodes2$celltype


# Map greenleaf colon cell type annotation ####
prj = 'greenleaf_colon'
archp_barcodes = rownames (archp@cellColData)[archp$project == prj]
archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '\\#'))[2])
#archp_barcodes = sapply (archp_barcodes, function(x) unlist (strsplit (x, '-'))[1])
archp_barcodes = data.frame (barcode = archp_barcodes, sample = archp$Sample[archp$project == prj])
archp_barcodes = split (archp_barcodes, archp_barcodes$sample)

ext_meta = meta[[6]]
ext_meta$barcode = sapply (ext_meta$Cell, function(x) unlist (strsplit (x, '\\#'))[2])
ext_meta$sample = ext_meta$Sample
ext_meta = split (ext_meta, ext_meta$sample)

match_sample = sapply (archp_barcodes, function(x) sapply (ext_meta, function(y) sum (x$barcode %in% y$barcode)))
names(archp_barcodes) = rownames(match_sample)[apply (match_sample, 2, which.max)]
archp_barcodes = archp_barcodes[!names(archp_barcodes) == '']
archp_barcodes2 = lapply (names(archp_barcodes), function(x) 
  {
  archp_barcodes[[x]]$celltype = paste0(ext_meta[[x]]$GrossPathology,'_',ext_meta[[x]]$CellType)[match(archp_barcodes[[x]]$barcode, ext_meta[[x]]$barcode)]
  archp_barcodes[[x]]$barcode = rownames(archp_barcodes[[x]])
  archp_barcodes[[x]]
  })
archp_barcodes2 = do.call (rbind, archp_barcodes2)
archp$celltype[match(archp_barcodes2$barcode, rownames(archp@cellColData))] = archp_barcodes2$celltype

archp$celltype_prj = paste0(archp$project, '_',archp$celltype)

projects = unique(archp$project)
#projects = c('yang_kidney','rawlins_fetal_lung')

archp = archp[!is.na(archp$celltype)]
