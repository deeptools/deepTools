computeMatrix reference-point -a 100 -b 100 -S test.bw -R test2.bed -o master.mat.gz -bs 1 -p 1
# unzip but keep original gz file.
gunzip -f -c master.mat.gz > master.mat

# test referencePoint center
computeMatrix reference-point -a 100 -b 100 --referencePoint center -S test.bw -R test2.bed -o master_center.mat.gz -bs 1
# unzip but keep original gz file.
gunzip master_center.mat.gz

# test referencePoint TES
computeMatrix reference-point -a 100 -b 100 --referencePoint center -S test.bw -R test2.bed -o master_TES.mat.gz -bs 1
# unzip but keep original gz file.
gunzip  master_center_TES.mat.gz

computeMatrix reference-point -R test2.bed -S test.bw  -b 100 -a 100 --outFileName master_nan_to_zero.mat.gz -bs 1 -p 1 --missingDataAsZero
gunzip -c  master_nan_to_zero.mat.gz > master_nan_to_zero.mat

computeMatrix scale-regions -a 100 -b 100 -m 100 -S test.bw -R test2.bed -o master_scale_reg.mat.gz -bs 1 -p 1
gunzip -c master_scale_reg.mat.gz > master_scale_reg.mat

plotHeatmap -m master.mat.gz --outFileName master.png
plotHeatmap -m master.mat.gz --outFileName master_relabeled.png --regionsLabel uno dos
plotHeatmap -m master_scale_reg.mat.gz --outFileName master_scale_reg.png

plotProfile -m master.mat.gz --outFileName profile_master.png --regionsLabel uno dos --plotType std
plotProfile -m master.mat.gz --outFileName profile_master_heatmap.png --plotType heatmap

# for tests with multiple bigwigs and multiple beds
computeMatrix reference-point -R group1.bed group2.bed -S test.bw  -b 100 -a 100 --outFileName master_multibed.mat.gz  -bs 1 -p 1
gunzip -c  master_multibed.mat.gz > master_multibed.mat

computeMatrix reference-point -R group1.bed group2.bed -S test.bw  -b 100 -a 500 --outFileName master_extend_beyond_chr_size.mat.gz -bs 1 -p 1
gunzip -c master_extend_beyond_chr_size.mat.gz > master_extend_beyond_chr_size.mat

computeMatrix reference-point -R group1.bed group2.bed -S test.bw test.bw test.bw test.bw -o master_multi.mat.gz -a 100 -b 100 -bs 1

plotHeatmap -m master_multi.mat.gz --perGroup --outFileName heatmap_master_multi_pergroup.png --samplesLabel file1 file2 file3 file4
plotHeatmap -m master_multi.mat.gz --colorList 'white,blue' 'white, red' --zMin 1 0 --zMax 4 5 -o heatmap_master_multi_color.png
plotHeatmap -m master_multi.mat.gz --colorMap Reds binary terrain --boxAroundHeatmaps no -o heatmap_master_multi_colormap_no_box.png
plotHeatmap -m large_matrix.mat.gz --interpolation bilinear --outFileName heatmap_master_interpolation_bilinear.png

plotProfile -m master.mat.gz --outFileName profile_master_overlap_lines.png --plotType overlapped_lines --yMin -1
plotProfile -m master_multi.mat.gz --outFileName profile_master_multi.png --numPlotsPerRow 2 --yMax 1.5
plotProfile -m master_multi.mat.gz --outFileName profile_master_multi_pergroup.png --perGroup --yMax 1.5

