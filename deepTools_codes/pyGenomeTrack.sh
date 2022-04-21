make_tracks_file --trackFiles /Users/jplab/Desktop/cut_tag_20210611/bigwig/FOXA1_B14.bw /Users/jplab/Desktop/cut_tag_20210611/bigwig/FOXA1_DMSO.bw /Users/jplab/Desktop/cut_tag_20210611/bigwig/FOXA1_TGF.bw /Users/jplab/Desktop/cut_tag_20210611/peak/lane_1/FOXA1_B14_pooled_peaks.narrowPeak -o FOXA1_cuttag_bw.ini



pyGenomeTracks --tracks ./FOXA1_cuttag_bw.ini --region chr14:75,936,001-76,004,694 --outFileName FOXA1_B14_bw_TGFB3.pdf