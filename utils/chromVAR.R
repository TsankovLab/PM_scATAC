archp = addBgdPeaks (archp, force= force)
  archp = addMotifAnnotations (ArchRProj = archp, 
      motifSet = "cisbp", 
      #motifSet = 'JASPAR2020',
      #name = "JASPAR2020_Motif",
      force=force)
  archp = addDeviationsMatrix (
    ArchRProj = archp, 
    peakAnnotation = "Motif",
    force = force
  )
  
  archp = saveArchRProject (ArchRProj = archp,  
      load = force)