<h1 align="center">T4</h1>
<h3 align="center">Calculates T4scores for each TTTT motifs in genome</h3>

<p align="left"> <img src="https://komarev.com/ghpvc/?username=13rajendra&label=Profile%20views&color=0e75b6&style=flat" alt="13rajendra" /> </p>

<img width="222" alt="Screen Shot 2024-02-14 at 5 12 47 PM" src="https://github.com/13rajendra/T4/assets/130776338/dd532a4f-81f7-4437-a85b-c0ddb2d1eaac">


This is a github package which does the following:

i) Find TTTT motifs across the genome using seqkit tool.

ii) Take 3' end signal of reads from Nascent-Seq and find the enrichment of signal over each of these TTTT motifs which we called as T4score.
          (_3'end signal of reads from nascent seq is already provided in ./output/_)

iii) Assign the maximum T4score that intersects our set of genes of interest.


