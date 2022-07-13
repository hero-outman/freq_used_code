color1='#f7dc6f'
color2='#F08080'
color3='#45B39D'
color4='#82E0AA'
color5='#A569BD'
color6='#2E86C1'

matrix=shFOXA1_ncFOXA1.TGF.groupByBed_ncshT.navieoverlap.center.3k.mat.gz 
matrix_basename=`basename $matrix .mat.gz`

# one bigwig, one plot
plotProfile \
    -m ${matrix} \
    --outFileName ${matrix_basename}.profile.png \
                   --plotType fill \
                   --colors $color1 $color2 $color3 $color4 $color5 $color6 \
                   --numPlotsPerRow 2 \
                   --regionsLabel 'peak intensity' \
                   --samplesLabel "ncFOXA1.B14" "shFOXA1.B14" "ncFOXA1.DMSO" "shFOXA1.DMSO" "ncFOXA1.TGF" "shFOXA1.TGF" \
                   --plotTitle 'peak intensity profile' \
                   --plotFileFormat png

# one bed, one plot
plotProfile \
    -m ${matrix} \
    --outFileName ${matrix_basename}.profile.overlay.png \
                   --plotType fill \
                   --colors $color1 $color2 $color3 $color4 $color5 $color6 \
                   --perGroup \
                   --regionsLabel 'ncFOXA1.B14 and shFOXA1.B14' 'ncFOXA1.DMSO and shFOXA1.DMSO' 'ncFOXA1.TGF-beta and shFOXA1.TGF-beta' \
                   --samplesLabel "ncFOXA1.B14" "shFOXA1.B14" "ncFOXA1.DMSO" "shFOXA1.DMSO" "ncFOXA1.TGF" "shFOXA1.TGF" \
                   --plotTitle 'peak intensity profile' \
                   --plotFileFormat png                   

plotProfile \
    -m ${matrix} \
    --outFileName ${matrix_basename}.profile.overlay.png \
                   --plotType fill \
                   --colors $color1 $color2 $color3 $color4 $color5 $color6 \
                   --perGroup \
                   --regionsLabel 'peak intensity' \
                   --samplesLabel "ncFOXA1.TGF" "shFOXA1.TGF" \
                   --plotTitle 'peak intensity profile' \
                   --plotFileFormat png         