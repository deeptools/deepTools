computeMatrix reference-point -a 100 -b 100 -S test.bw -R test2.bed -o master.mat.gz -bs 1
computeMatrix scale-regions -a 100 -b 100 -m 100 -S test.bw -R test2.bed -o master_scale_reg.mat.gz -bs 1
plotHeatmap -m master.mat.gz --outFileName master.svg
plotHeatmap -m master.mat.gz --outFileName master_relabeled.svg --regionsLabel uno,dos
plotHeatmap -m master_scale_reg.mat.gz --outFileName master_scale_reg.svg
plotProfile -m master.mat.gz -o profile_master.svg --regionsLabel uno,dos --plotType std
plotProfile -m master.mat.gz -o profile_master_heatmap.svg --plotType heatmap
