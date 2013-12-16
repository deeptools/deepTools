* [FAQ](#FAQ)

General hints about optimal deepTools usage
==========================================

* when you're playing around with the tools to see what kinds of results they will produce: limit the operation to one chromosome only to __save computation time__! ("advanced output options" --> "Region of the genome to limit the operation to")



Frequently asked questions
===========================

## How can I determine basic parameters of a BAM file, such as the number of reads, read length, duplication rate and average DNA fragment length?
Eventhough [MACS][] is meant to do peak calling for you, it also outputs a number of useful information such as those listed above.
Simply run MACS on the BAM file that you would like to gain the information for and check the .xls file from the MACS output. It will list:
* tag length = read length
* duplication rate
* number of tags = number of reads
* d = distance = average DNA fragment size

##
