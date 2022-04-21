import glob, os
from glbase3 import *

print(os.getcwd())
os.chdir('/Users/jplab/Desktop/cut_tag_20210611/filter_region/6-24')

config.draw_mode = "png"

hg38_path = '/Users/jplab/Desktop/cut_tag_20210611/filter_region/6-21/hg38_ensembl_v95_enst.glb'
hg38 = glload(hg38_path)

list_of_peaks = sorted(list(glob.glob("/Users/jplab/Desktop/cut_tag_20210611/filter_region/6-23/B14_FOXA1.ATACUniq.K27ac.bed")))

for filename in list_of_peaks:
    print('filename: '+filename)
    gl = genelist(filename, format=format.bed)
    gl.name = os.path.split(filename)[1].replace(".bed", "")
    print('glname: '+ gl.name)
    print('gl: '+ gl)

    # reverse the order of the list, in place
    gl.reverse()

    print('reversed gl: ', gl)

    # image_filename = gl.name +'.png'
    # ann = hg38.annotate(genelist=gl, image_filename=image_filename, distance=1000, closest_only=False,window=500)
    # # ann = ann.removeDuplicates("ensg")
    # ann.saveTSV("./%s-ann1000kb.tsv" % (gl.name), key_order=['ensg', 'enst'])
    # ann.save("./%s-ann1000kb.glb" % (gl.name))

    # image_filename = gl.name +'.5k.png'
    # ann = hg38.annotate(genelist=gl, image_filename=image_filename, distance=5000, closest_only=False,window=500)
    # # ann = ann.removeDuplicates("ensg")
    # ann.saveTSV("./%s-ann5000kb.tsv" % (gl.name), key_order=['ensg', 'enst'])
    # ann.save("./%s-ann5000kb.glb" % (gl.name))

    # ann = hg38.annotate(genelist=gl, image_filename=image_filename, distance=5000, closest_only=False,window=500)
    # ann = ann.removeDuplicates("ensg")
    # ann.saveTSV("./%s-ann5000kb-nonamedupes.tsv" % (gl.name), key_order=['ensg', 'enst'])
    # ann.save("./%s-ann5000kb-nonamedupes.glb" % (gl.name))

    # ann = hg38.annotate(genelist=gl, image_filename=image_filename, distance=5000, closest_only=True,window=500)
    # ann = ann.removeDuplicates("ensg")
    # ann.saveTSV("./%s-ann5000kb-nonamedupes.closest.tsv" % (gl.name), key_order=['ensg', 'enst'])
    # ann.save("./%s-ann5000kb-nonamedupes.closest.glb" % (gl.name))

    # ann = hg38.annotate(genelist=gl, image_filename=image_filename, distance=1000, closest_only=False,window=500)
    # ann = ann.removeDuplicates("ensg")
    # ann.saveTSV("./%s-ann1000kb-nonamedupes.tsv" % (gl.name), key_order=['ensg', 'enst'])
    # ann.save("./%s-ann1000kb-nonamedupes.glb" % (gl.name))
