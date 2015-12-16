computeMatrix reference-point -a 100 -b 100 -S test.bw -R test2.bed -o master.mat.gz -bs 1
computeMatrix scale-regions -a 100 -b 100 -m 100 -S test.bw -R test2.bed -o master_scale_reg.mat.gz -bs 1
plotHeatmap -m master.mat.gz --outFileName master.svg
plotHeatmap -m master.mat.gz --outFileName master_relabeled.svg --regionsLabel uno,dos
plotHeatmap -m master_scale_reg.mat.gz --outFileName master_scale_reg.svg
plotProfile -m master.mat.gz -o profile_master.svg --regionsLabel uno,dos --plotType std
plotProfile -m master.mat.gz -o profile_master_heatmap.svg --plotType heatmap
plotProfile -m master.mat.gz -o profile_master_overlap_lines.svg --plotType overlapped_lines --yMin -1

# for tests with multiple bigwigs and multiple beds
computeMatrix reference-point -R group1.bed group2.bed -S test.bw test.bw test.bw test.bw -o master_multi.mat.gz -a 100 -b 100 -bs 1
plotHeatmap -m master_multi.mat.gz -o heatmap_master_multi_pergroup.svg --perGroup
plotProfile -m master_multi.mat.gz -o profile_master_multi.svg --numPlotsPerRow 2 --yMax 1.5
plotProfile -m master_multi.mat.gz -o profile_master_multi_pergroup.svg --perGroup --yMax 1.5
