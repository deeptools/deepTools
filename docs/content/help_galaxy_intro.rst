Using deepTools within Galaxy    
================================

We have a publicly available deepTools installation embedded within the Galaxy framework: `deeptools.ie-freiburg.mpg.de <http://deeptools.ie-freiburg.mpg.de/>`__

This server also contains some additional tools that will enable users to analyse and visualize data from high-throughput sequencing experiments 

The information on this page are meant to: 

* **provide a brief introduction into the Galaxy framework in general:** `Galaxy <http://galaxyproject.org/>`__ is a tremendously useful platform developed by the Galaxy Team at Penn State and the Emory University. This platform is meant to offer access to a large variety of bioinformatics tools that can be used without computer programming experiences. That means, that the basic features of Galaxy will apply to every tool, i.e. every tool provided within a Galaxy framework will look very similar and will follow the concepts of Galaxy.

* **introduce you to the specific tools that are available in our web server dedicated to NGS data analysis**

  If you do not know the difference between a BAM and a BED file, that's fine, you can read up on them in our :doc:`help_glossary`.

For more specific help, check our :doc:`help_faq_galaxy` and the :doc:`example_step_by_step`.

**Table of contents**

.. contents:: 
    :local:

Basic features of Galaxy
-------------------------

The Galaxy team develops the platform, but since it is impossible to
meet all bioinformatics needs (that can range from evolutionary analysis
to data from mass spectrometry to high-throughput DNA sequencing (and
way beyond)) with one single web server, many institutes have installed
their own versions of Galaxy tuned to their specific needs. Our
`deepTools Galaxy <http://deeptools.ie-freiburg.mpg.de/>`__ is such a
specialized server dedicated to the analysis of high-throughput DNA
sequencing data. The overall makeup of this web server, however, is the
same as for any other Galaxy installation, so if you've used Galaxy
before, you will learn to use deepTools in no time!

The start site
^^^^^^^^^^^^^^^

Here is a screenshot of what the start site will approximately look like: 

.. image:: ../images/Gal_startsite.png

The start site contains 4 main features:

Top menu
~~~~~~~~~~~~~~

This will lead you to other sections of Galaxy, away from the actual analysis part.
Registered users, for example, can create and access **workflows**.
You can also find **content shared with you** by other users such as sample data sets, pages and workflows.

Tool panel
~~~~~~~~~~~~~~

*"What can be done"*: This panel lists all the *tools* installed in this Galaxy instance.

Main frame
~~~~~~~~~~~~~~

*"What am I doing?"*: the center frame is your main **working space** where input will be required from you once you've selected a specific tool.
In addition, you will always find **general information about the tool** here.

History panel
~~~~~~~~~~~~~~

*"What did I do?"*: here you can find all *files* that one produces or uploads

The history is like a **log book**: everything you ever did is recorded here (unless you deleted things permanently)

- histories can be shared with other users
- histories can also be downloaded
- for each file that was produced, you will find all kinds of useful  information such as the tool that was used to create the file, the tool's parameters etc.

For those visual learners, here's an annotated screenshot:

.. image:: ../images/Gal_startsite_with_comments.png

In the default state of the tool panel you see the **tool categories**,
e.g. "Get Data". If you click on them, you will see the **individual
tools** belonging to each category, e.g. "Upload File from your
computer", "UCSC Main table browser" and "Biomart central server" in
case you clicked on "Get Data". To use a tool such as "Upload File from
your computer", just click on it.

The **tool *search* panel** is extremely useful as it allows you to
enter a key word (e.g. "bam") that will lead to all the tools mentioning
the key word in the tool name.

Once you've uploaded any kind of data, you will find the history on the
right hand side filling up with green tiles. Each tile corresponds to
one data set that you either uploaded or created. The data sets can be
images, raw sequencing files, text files, tables - virtually anything.
The content of a data set *cannot* be modified - every time you want to
change something *within* a data file (e.g. you would like to sort the
values or add a line or cut a column), you will have to use a Galaxy
tool that will lead to a *new* data set being produced. This behaviour
is often confusing for Galaxy novices (as histories tend to accumulate
data sets very quickly), but it is necessary to enforce the strict
policy of documenting *every modification* to a given data set.
Eventhough your history might be full of data sets with strange names,
you will always be able to track back the source and evolution of each
file. Also, every data set can be downloaded to your computer.

Have a look at the following screenshot to get a feeling for how many
information Galaxy keeps for you (which makes it very feasible to
reproduce any analysis):

.. image:: ../images/Gal_screenshot_dataSet.png

Each data set can have 4 different states that are intuitively
color-coded:

.. image:: ../images/Gal_screenshot_dataSetStates.png

.. tip:: If you encounter a failure after you've run a tool, please follow the steps (in this order):

         1. click on the center button on the lower left corner of the failed data set ``(i)``: did you chose the **correct data files**?
         2. if you're sure that you chose the correct files, hit the ``re-run button`` (blue arrow in the lower left corner) - check again whether your files had the **correct file format**. If you suspect that the format might be incorrectly assigned (e.g. a file that should be a BED file is labelled as a tabular file), click the ``edit button`` (the pencil) of the input data file - there you can change the corresponding attributes
         3. if you've checked your input data and the error is persisting, click on the ``green bug`` (lower left corner of the failed data set) and send the **bug report** to us. You do not need to indicate a valid email-address unless you would like us to get in touch with you once the issue is solved.

Data import into Galaxy
-------------------------

There are three main ways to populate your Galaxy history with data
files:

.. #############################################################
   MAKE INTERNAL LINKS HERE ###
   OR MAKE A NEW PAGE
   ############################################################

1. Data upload from your computer
2. Import a shared data set from the Galaxy data library
3. Download annotation data from public servers

an additional option is to Copy data sets between histories


Upload files from your computer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The data upload of files **smaller than 2 GB** that lie on your computer is fairly straight-forward: click on the category "Get data" and choose the tool "Upload file".
Then select the file via the "Browse" button.

.. image:: ../images/Gal_DataUpload.png

For files **greater than 2GB**, there's the option to upload via an FTP server. If your data is available via an URL that links to an FTP server, you can simply
paste the URL in the empty text box.

If you do not have access to an FTP server, you can directly upload to our Galaxy's FTP.

 * first register with deeptools.ie-freiburg.mpg.de (via “User” ⟶ “register”; registration requires an email address and is free of charge)
 * You will also need an FTP client, e.g. `filezilla <https://urldefense.proofpoint.com/v2/url?u=https-3A__filezilla-2Dproject.org_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=qIjL9RRxwt_RObaavha0-3PJavlW5JAAePP8g6_zRFM&e=>`__.
 * Then login to the **FTP client** using your **deepTools Galaxy user name and password** (host: deeptools.ie-freiburg.mpg.de). Down below you see a screenshot of what that looks like with filezilla.
 * Copy the file you wish to upload to the remote site (in filezilla, you can simply drag the file to the window on the right hand side)
 * Go back to `deepTools Galaxy <https://urldefense.proofpoint.com/v2/url?u=http-3A__deeptools.ie-2Dfreiburg.mpg.de_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=1xpNM-7I4Z6-ZIJErjnO726mjBKdGH92RCWOc5kGh-U&e=>`__
 * Click on the tool "Upload file" (⟶ "Files uploaded via FTP") - here, the files you just copied over via filezilla should appear. Select the files you want and hit “execute”. They will be moved from the FTP server to your history (i.e. they will be deleted from the FTP once the upload was successful).

.. image:: ../images/Gal_filezilla.png

Import data sets from the Galaxy data library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you would like to play around with sample data, you can import files
that we have saved within the general data storage of the deepTools
Galaxy server. Everyone can import them into his or her own history,
they will not contribute to the user's disk quota.

You can reach the data library via "Shared Data" in the top menu, then
select "Data Libraries".

Within the Data Library you will find a folder called "Sample Data" that
contains data that we downloaded from the `Roadmap project <https://urldefense.proofpoint.com/v2/url?u=http-3A__www.roadmapepigenomics.org_data&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=th-6vqsGlAXEh96RIzGXdL-u2ypvcD6g-BA86le-Y5A&e=>`__ and
`UCSC <https://urldefense.proofpoint.com/v2/url?u=http-3A__genome.ucsc.edu_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=-tpAnqa6rqhffRRpItAYHCwcIb0KVxgd4jW667tchPk&e=>`__. More precisely, we downloaded the
:ref:`fastq` files and mapped the reads to the human reference genome
(version hg19) to obtain the :ref:`BAM <bam>` files you see. In addition, you will
find signal tracks of DNase-seq data from UCSC, bigWig files with GC
content for flies and mice and some annotation files.

Download annotation files from public data bases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many cases you will want to query your sequencing data results for
known genome annotation, such as genes, exons, transcription start sites
etc. These information can be obtained via the two main sources of
genome annotation, `UCSC <https://urldefense.proofpoint.com/v2/url?u=http-3A__genome.ucsc.edu_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=-tpAnqa6rqhffRRpItAYHCwcIb0KVxgd4jW667tchPk&e=>`__ and
`BioMart <https://urldefense.proofpoint.com/v2/url?u=http-3A__www.biomart.org_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=Et46CytirsKJYRV1jwPA3FSVUzJBAHLWYJUcOuHjBYQ&e=>`__. Please note that UCSC and BioMart
will cater to different ways of genome annotation, i.e. genes defined in
UCSC might not correspond to the same regions in a gene file downloaded
from BioMart. (For a brief overview over the issues of genome
annotation, you can check out
`Wikipedia <https://urldefense.proofpoint.com/v2/url?u=http-3A__en.wikipedia.org_wiki_Genome-5Fproject&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=40yftLP0yJfyWIQ056g47LjDbHES8e6UnOaMO9dPhCo&e=>`__, if you'd
always wanted to know much more about those issues,
`this <http://www.ncbi.nlm.nih.gov/pubmed/22510764>`__ might be a good
start.)

You can access the data stored at UCSC or BioMart conveniently through
our Galaxy instance which will import the resulting files into your
history. Just go to **"Get data"** ⟶ "UCSC" or "BioMart".

The majority of annotation files will probably be in BED format,
however, you can also find other data sets. UCSC, for example, offers a
wide range of data that you can browse via the "group" and "track" menus
(for example, you could download the GC content of the genome as a
signal file from UCSC via the "group" menu ("Mapping and Sequencing
Tracks"). Note, however, that the download through this interface is
limited to 100,000 lines per file which might not be sufficient for some
mammalian data sets).

Here's a screenshot from downloading a BED-file of all RefSeq genes
defined for the human genome (version hg19):

And here's how you would do it for the BioMart approach:

Per default, **BioMart will not output a BED file** like UCSC does. It
is therefore important that you make sure you get all the information
you need (most likely: chromosome, gene start, gene end, ID, strand) via
the "Attributes" section. You can click on the "Results" button at any
time to check the format of the table that will be sent to Galaxy (Note
that the strand information will be decoded as 1 for "forward" or "plus"
strand and -1 for "reverse" or "minus" strand.)

    Be aware, that BED files from UCSC will have chromosomes labelled
    with “chr” while ENSEMBL usually returns just the number – this
    might lead to incompatibilities, i.e. when working with annotations
    from UCSC and ENSEMBL, you need to make sure to use the same naming!

 #### Copy data sets between histories In case you have registered with
deepTools Galaxy you can have more than one history. In order to
minimize the disk space you're occupying we strongly suggest to **copy**
data sets between histories when you're using the same data set in
different histories. This can easily be done via the History panel's
option button ⟶ "Copy dataset". In the main frame, you should now be
able to select the history you would like to copy from on the left hand
side and the target history on the right hand side.

`Back to the deepTools Galaxy <https://urldefense.proofpoint.com/v2/url?u=http-3A__deeptools.ie-2Dfreiburg.mpg.de_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=1xpNM-7I4Z6-ZIJErjnO726mjBKdGH92RCWOc5kGh-U&e=>`__

Which tools can I find in the deepTools Galaxy?
-----------------------------------------------

As mentioned above, each Galaxy installation can be tuned to the
individual interests. Our goal is to provide a Galaxy that enables you
to **quality check, process and normalize and subsequently visualize
your data obtained by high-throughput DNA sequencing**.

We provide the following kinds of tools:

1. :ref:`deepTools - NGS data handling <help_galaxy_intro_deepTools>`
2. :ref:`peak calling (ChIP-seq specific) <help_galaxy_intro_peaks>`
3. :ref:`operating on genomic intervals <help_galaxy_intro_BED>`
4. :ref:`working with text files and tables <help_galaxy_intro_textfiles>`

.. _help_galaxy_intro_deepTools:

deepTools
^^^^^^^^^

The most important category is called **"deepTools"** that contains 8
major tools (for information on the data formats, see our
:doc:`help_glossary`)

We have compiled several sources of detailed information specifically
about the usage of deepTools:

1. General overview of :doc:`how we use deep Tools <example_usage>`
2. Each individual tool is described in more detail on separate pages :doc:`here <list_of_tools>`
3. For each tool, you will find specific explanations within the
   `deepTools Galaxy <https://urldefense.proofpoint.com/v2/url?u=http-3A__deeptools.ie-2Dfreiburg.mpg.de_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=1xpNM-7I4Z6-ZIJErjnO726mjBKdGH92RCWOc5kGh-U&e=>`__ main
   frame, too.
.. 4. the `example workflows <Example-workflows>`__ might help to get a
   feeling for the kinds of analyses than can be done with `deepTools
   Galaxy <https://urldefense.proofpoint.com/v2/url?u=http-3A__deeptools.ie-2Dfreiburg.mpg.de_&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=1xpNM-7I4Z6-ZIJErjnO726mjBKdGH92RCWOc5kGh-U&e=>`__

.. _help_galaxy_intro_peaks:

Peak calling
^^^^^^^^^^^^^^

In ChIP-seq analysis, peak calling algorithms are essential downstream
analysis tools to identify regions of significant enrichments (i.e.
where the ChIP sample contained significantly more sequenced reads than
the input control sample). By now, there must be close to 100 programs
out there (see `Wilbanks et
al. <https://urldefense.proofpoint.com/v2/url?u=http-3A__www.plosone.org_article_info-253Adoi-252F10.1371-252Fjournal.pone.0011471&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=lhLQ7qst_E5ZweBT_PdS_mJIE9biseGu2DTBPk2papM&e=>`__
for a comparison of peak calling programs).

In contrast to deepTools that were developed for handling and generating
*continuous* genome-wide profiles, peak calling will result in a *list
of genomic regions*. Have a look at the screenshot to understand the
difference.

We have included the peak callers
`MACS <http://www.ncbi.nlm.nih.gov/pubmed/22936215>`__ and
`SICER <https://urldefense.proofpoint.com/v2/url?u=http-3A__bioinformatics.oxfordjournals.org_content_25_15_1952.full&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=4ZEmdk9_IT-qF0ZDdKMF6Z-vWNUrYB3r76ucGWLaCYo&e=>`__
within our Galaxy instance with
`MACS <http://www.ncbi.nlm.nih.gov/pubmed/22936215>`__ being the most
popular peak calling algorithm for the identification of localized
transcription factor binding sites while
`SICER <https://urldefense.proofpoint.com/v2/url?u=http-3A__bioinformatics.oxfordjournals.org_content_25_15_1952.full&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=4ZEmdk9_IT-qF0ZDdKMF6Z-vWNUrYB3r76ucGWLaCYo&e=>`__
was developed for diffuse ChIP-seq signals. Note that MACS version 1.14
is quite different from MACS version 2.

.. _help_galaxy_intro_BED:

Working with genomic intervals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Galaxy has 2 file formats to store lists of genomic regions:

-  INTERVAL

   -  tab-separated
   -  requirements:

      1. Column: chromosome
      2. Column: start position
      3. Column: end position

   -  all other columns can contain any value or character

-  BED

   -  very similar to INTERVAL, but stricter when it comes to what is
      expected to be kept in which column:

      -  

         1. to 3. Column: same as interval

      -  Column 4: name
      -  Column 5: score
      -  Column 6: strand

In case you would like to work with several lists of genomic regions,
e.g. generate a new list of regions that are found in two different
files etc., there are two categories of tools dedicated to performing
these tasks:

 * Operate on genomic intervals
 * BEDtools

Each tool's function is explained within Galaxy. Do browse those tools
as they will give you a very good glimpse of the scope of possible
analyses!

.. _help_galaxy_intro_textfiles:

Working with text files and tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to deepTools that
were specifically developed for the handling of NGS data, we have
incorporated several standard Galaxy tools that enable you to manipulate
tab-separated files such as gene lists, peak lists, data matrices etc.

There are 3 main categories:

-  **Text manipulation**

   -  unlike Excel where you can easily interact with your text and
      tables via the mouse, data manipulations within Galaxy are
      strictly based on commands. If you feel like you would like to do
      something to certain *columns* of a data set, go through the tools
      of this category
   -  e.g. adding columns, cutting columns, pasting two files side by
      side, selecting random lines etc.
   -  a very useful tool of this category is called *Trim* - if you need
      to remove some characters from a column, this tool's for you! (for
      example, sometimes you need to adjust the chromosome naming
      between two files from different source - using *Trim*, you can
      remove the "chr" infront of the chromosome name)

-  **Filter and Sort**

   -  in addition to the common sorting and filtering, there's the very
      useful tool to *select lines that match an expression* (for
      example, using the expression *c1=='chrM'* will select all rows
      from a BED file with regions located on the mitochondrial
      chromosome)

-  **Join, Subtract, Group**
-  this category is very useful if you have several data sets that you
   would like to work with, e.g. by comparing them

Workflows
--------------------

Workflows are Galaxy's equivalent of protocols.
This is a very useful feature as it allows users to *share
their protocols and bioinformatic analyses* in a very easy and
transparent way.
This is the graphical representation of a Galaxy
workflow that can easily be modified via drag'n'drop within the
workflows manual (you must be registered with deepTools Galaxy to be
able to generate your own workflows or edit published ones).

Where to get help?
--------------------

Please check our `deepTools Galaxy FAQs <Galaxy-related-FAQs>`__

-  `general Galaxy help <https://urldefense.proofpoint.com/v2/url?u=http-3A__wiki.galaxyproject.org_Learn&d=BQIGaQ&c=lb62iw4YL4RFalcE2hQUQealT9-RXrryqt9KZX2qu2s&r=YPs4H2QfvX0QdeqqpLIqoKZMYe9vwL5KkadTIhRrkBU&m=V0hrMSIcFCpE37KzRB4Nzvnu1qyvX8PcXgnmi5X4OxU&s=8jffwhatG2bweZYXURXCIXFA05QO7BfwfAFQdwE3azc&e=>`__
-  specific help with deepTools Galaxy: deeptools@googlegroups.com
-  if you encounter a failing data set (marked in red), please send a
   bug report via Galaxy and we will get in touch

