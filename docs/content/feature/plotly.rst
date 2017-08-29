Plotly
======

`Plotly <http://plot.ly>`__ is a tool for interactive visualization of datasets that deepTools has supported since version 2.6. The output of this "image format" is a web page that you can load in your browser of choice. Generally the resulting webpages will look essentially like the default images, but you can do things such as mouse-over an area to get the sample label or value and zoom in/pan around images. These web pages can grow to be **very** large, particularly for `plotHeatmap`. Further, not all deepTools options will be supported. In particular, in `plotHeatmap`, only a single colorbar is supported. These limitations are due to limitations in plotly itself and may be removed as plotly itself matures.

Note that plotly output is saved locally and remote servers are not used. This is due to most users likely wishing to keep their data private. If you would like to edit the resulting pages and export them to plot.ly then click on the link in the bottom-right corner of the plotly output.
