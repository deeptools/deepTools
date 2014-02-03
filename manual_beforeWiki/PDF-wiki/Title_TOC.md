# deepTools: a flexible platform for exploring deep-sequencing data
---------------------------
----------------------------
## MANUAL

1. **Why we built deepTools**
2. **How we use deepTools**
3. **What deepTools can do**
4. **Tool details**
	* Quality controls of aligned reads
		* _bamCorrelate_
		* _computeGCbias_
		* _bamFingerprint_
	* Normalization and bigWig generation
		* _correctGCbias_ 
		* _bamCoverage_
		* _bamCompare_ 
	* Visualization: heatmaps and summary plots
5. **Glossary: Abbreviations and file formats**

![figI](https://raw.github.com/fidelram/deepTools/master/examples/collage_wout_header.png)


| Fidel Ramírez, Friederike Dündar, Sarah Diehl, Björn A. Grüning, Thomas Manke |
|:-------------------:|
| _Bioinformatics Group, Max-Planck-Institute of Immunobiology and Epigenetics & Department of Computer Science, University of Freiburg_ |

Web server (incl. sample data): [deepTools.ie-freiburg.mpg.de](http://deepTools.ie-freiburg.mpg.de)
Code: [github.com/fidelram/deepTools](https://github.com/fidelram/deepTools/)