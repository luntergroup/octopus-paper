# Reviewer #1:

Remarks to the Author:
In their manuscript, Cooke, Wedge, and Lunter describe Octopus, their haplotype-based software tool for variant calling in next-generation sequencing (NGS) data. They compare their method to several widely used tools across multiple experimental designs, including germline calling, de novo mutation detection in trios, and somatic mutation detection with paired/unpaired samples. Octopus consistently performs comparably to or better than these tools in the scenarios presented. 

In addition to reviewing the manuscript, I tested its code and sample data through the Code Ocean container platform. I found the files and documentation provided with the package to be thorough, and was able to run the program without issue. 

Overall, I found the analyses presented to be suitable for this type of manuscript. Major comments are as follows:

1.) The authors write in the abstract that "Octopus accurately characterizes germline and somatic variation in tumors, both with and without a paired normal sample." This last portion of this statement already raises some alarm, as numerous have reported on, and attempted to address, the significant challenge of calling somatic mutations in tumor samples that lack matched normal samples. In the authors' hands, the accuracy of their tool (F-measure for the somatic test in Supplementary Table 3) was 0.7199 for the synthetic skin tumor and 0.4094 for the synthetic breast tumor. While Octopus out-performs another tool (Pisces) in this feat, this is still somewhat poor performance overall. My concern is that this might mislead readers into believing that Octopus magically solves the significant problem of somatic-calling in tumor-only samples. This concern is compounded by the decision to compute "overall" accuracy that includes germline variants as part of the same analysis. I find
that fundamentally problematic, since identifying the overall variant complement in tumor-only samples is not a significant challenge.

<span style="color:red">
Just re-word abstract?
</span>

2.) The authors highlight microinversion detection as a feature of Octopus. This class of variation is not usually considered by most variant callers, so it would represent a true advance in capability. However, according to my reading of the online methods (page 10), "Complex variation such as microinversions are found by inspecting bubbles" that are produced by local reassembly. To me, this means that a human must review local microassemblies and manually call microinversions. If so, then Octopus does not detect these variations and it is misleading to claim such (it's akin to saying that IGV calls SNPs, because a human can use the viewer to look at aligned BAM files and call variants). 

<span style="color:red">
Reviewer is misunderstanding. Octopus outputs microinversions calls explicetly. e.g.:

```shell
2	134937085	.	CAGAGACTTTATAGAACAAGAC	GTCTTGTTCTATAAAGTCTCTG	2783.56	PASS	AC=2;AN=2;DP=47;MQ=60;MQ0=0;NS=1;RFQUAL=26.6508	GT:GQ:DP:MQ:PS:PQ:FT	1|1:146:47:60:134937015:99:PASS
```
</span>

3.) On page 3, the authors suggest that "Octopus is more robust to noise than other methods" based on its performance on two public datasets generated on the Illumina HiSeq X Ten instrument. A different error profile is only one of the characteristics that separate X Ten data from that of other sequencing instruments. In my opinion, the authors have not provided sufficient evidence to make this suggestion about their caller being more robust to noise. 

<span style="color:red">
Fair enough... remove sentence?
</span>

4.) I have a few issues with the de novo mutation calling analysis described on page 4. First, the authors compare Octopus running in trio mode -- which presumably is designed specifically for detecting de novo mutations with confidence -- to a number of other tools in joint-calling mode (i.e. calling the trio like a cohort of 3 samples). This seems like an apples-to-oranges comparison, and I'd further point out that there are tools like DeNovoGear that are specifically designed for de novo mutation detection in trios, that might be better suited for this comparison. Secondly, the authors base their claims on a single WGS trio (confusingly named WGS500). This is not nearly enough data to assess performance in de novo mutation calling. If the authors wish to highlight their tool's abilities, they need to provide a more substantial dataset and compare it with the right tools.

<span style="color:red">
Most published de novo studies use custom pipelines after GATK joint calling. Finding good validation data difficult. Worth doing some sequencing?
</span>

5.) On page 6, the authors claim that "Overall, Octopus had substantially higher F-measure on both tests than all other methods (Supplementary Table 2)." I'll point out that Supplementary Table 2 includes the results of the downsampling experiments, in some of which other callers returned higher F-measures than Octopus. Thus, I feel this is a bit of an overclaim. In the two synthetic tumor data sets, Octopus returned a higher F-measure than other tools. In MOST but not all downsampling experiments, it did the same.

<span style="color:red">
We get higher F-Measure on all tests when using RFQUAL threshold of 7 (although this decreases main test F-Measure slightly).
</span>

6.) With regard to phasing, which is a promising capability, Octopus only manages to phase 22% of variants in the paired synthetic skin tumor test. Because it and any other tool will be limited by such factors as read length, insert size, and the density of variation, it would be very useful if the authors explored phasing performance in more details. For example, I would like to know the average and distribution of distances to nearby variants to which Octopus succesfully phased a somatic variant. 

<span style="color:red">
Should be fairly straightforward...
</span>

7.) The authors remark that they "improved upon" the synthetic datasets generated by the ICGC-TCGA DREAM challenge. I would like to know the basis for the improvements and the justification for creating their own synthetic datasets rather than using an already-generated, well-studied one. The DREAM challenge datasets were configured for precisely the purpose of comparing variant caller performance, and a number of skilled teams from around the world applied their tools to the challenge of somatic mutation detection. The underlying data, truth sets, and top performer results are all publicly available, so I wonder why the authors did not provide a head-to-head comparison that benchmarks their tool on the same dataset.

<span style="color:red">
We explain at length why existing datasets are not ideal (and specifically mention DREAM data).
</span>

Minor points are as follows:
A.) In the Results section, paragraph 1, the authors write that "variants from existing VCF files may also be considered." Considered how? As priors, or as target variants to be characterized with the haplotype assembly? This should be clarified.
B.) Same paragraph, there is an agreement error: "Calls are made once there is sufficient confident."



# Reviewer #2:

Remarks to the Author:
Octopus is an open source, well documented haplotype-based variant caller that handles a wide variety of use cases. It especially tackles difficult problems like de novo calling, small indels and microdeletions in germline data, and somatic low frequency tumor-only calling and phasing. The authors do an nice job of validating across multiple sample types in the paper, covering all of the major tools and validation sets currently available. Thank you for making Octopus freely available and formally writing up the methods and validations.

I've previously used Octopus in comparison on low frequency UMI tagged tumor only data (https://github.com/bcbio/bcbio_validations/tree/master/somatic-lowfreq#vardict-156-octopus-051b) so have practical experience running Octopus and interacting with the authors. They've been responsive to suggestions and provided improvements and parameter tweaking suggestions. Below are my comments to help improve the manuscript, which center around:

- Some ideas for contacting tool authors for specific comparisons, where we might be able to better reflect other tool capabilities.
- Thoughts on data availability and additional tests for low frequency somatic variants.
- Request for rough timing estimates to place Octopus runtime in the context of other tools.

Validation suggestions:

- Strelka2 performance on 10X input data: it would be worthwhile getting in touch with Chris Saunders to see if the Illumina team has insight or comments for the paper. We looked in depth at 10X data for GATK4, unfortunately not Strelka2, but found that a combination of adapter trimming and filters helped improve calling:

https://github.com/bcbio/bcbio_validations/tree/master/gatk4#na24385-10x-data-on-grch37

<span style="color:red">
Chris has contacted us abd stated his only concern was with the de novo tests, which we will address.
</span>

- Similarly, it would be worth approaching the DeepVariant and Strelka2 development teams about the large number of de novo false positives with those callers. They might have suggestions for filtering or calling beyond what you've tried and be a more accurate representation for the paper.

<span style="color:red">
Andy Carroll from Mark DePristo's DeepVariant team has been in touch with us and acknolwedged, but was not able to offer suggestions for the de novo analysis. He did request that we update our results for the latest DeepVariant version which we have now done.
</span>

- The comments about phased indels and SNPs haplotypes in truth sets like GiaB are useful and would be great to relay back to that community to see if they can improve the underlying truth sets. It would be worth a proposal for how the variant community can better characterize these divergent representations.

<span style="color:red">
Is this a request?
</span>

- For the somatic mutation synthetic tumor, what was the cause of the incorrect spike in mutations from BAMsurgeon? I understood the cell line differences but not the underlying cause of these.

<span style="color:red">
Will look into this.
</span>

- Thank you for plotting recall versus variant allele frequency, this is the key component of somatic variant calling as we increasingly want to differentiate low frequency variants contributing to subclones. Have you done comparisons with higher depth samples beyond 60x and balanced tumor and normal sequencing depths? While downsampling is of interest, we see less of this with NovaSeqs now available, and more high depth panels where the goal is to identify <1% variants. It would be useful to understand how Octopus does in these cases.

<span style="color:red">
Should we re-do with higher depths?
</span>

- Are you able to share any use of Octopus on identifying phased quasispecies mutations in viral data? The text currently indicates it would be straightforward, although I've found that not to be true when evaluating on HIV data. I explored this a few years ago with several different low frequency callers available at the time (VarDict, LoFreq, FreeBayes) and had difficulty producing useful results, especially with any level of phasing (https://github.com/hbc/li_hiv_call3). I'd love to see any benchmark or usage suggestions from work you've done so far, as phased quasispecies results are important for HIV drug resistance predictions.

<span style="color:red">
Should we add bacteria results?
</span>

- In my validations, MuTect2 has improved greatly between the 4.0.0.0 release you evaluated and the current version. I know there are always new versions and validation is a never ending job, but it would be worth contacting the Broad MuTect2 team to see if they think your current results match their expectations.

<span style="color:red">
Updated results for newer Mutect2.
</span>

Paper suggestions:

- Please make the synthetic tumor data available. As you mention, we have a lack of good somatic datasets, and the community could provide comparisons of Octopus with the latest versions of MuTect2 and other tools as a complement to the work you've done.

<span style="color:red">
Need to discuss this.
</span>

- Could you provide rough reporting of runtimes for Octopus versus other methods on some major use cases (germline, trio, low frequency somatic)? An exhaustive comparison isn't necessary but in my experience there is a time tradeoff for using Octopus on low frequency cases versus other sensitive callers like VarDict, and it would be useful to quantify that in the paper or documentation.

<span style="color:red">
:(
</span>

- For the precision/recall curves, would it be possible to better distinguish the callers via different points? My color differentiation genetics is not great and it's hard to tell which caller is which other than the easily distinguishable Octopus line.

<span style="color:red">
Sure
</span>

Brad Chapman



# Reviewer #3:

Remarks to the Author:
This paper describes a new method for variant calling, implemented as a software package called Octopus, that the authors claim offers very good performance in a variety of scenarios. The authors do a nice job explaining the background, especially the issues for accurate germline variant calling and they are thorough in the sense that Octopus is compared to many established methods, including GATK, Strelka, and FreeBayes. However, the scientific advance provided by this work is limited.

There are a number of criticisms. 

First, some results are overstated, for example Octopus having "the highest F-measure on all tests other than...HG005" in Fig 2. Differences with several callers are small and, in some cases, virtually negligible. They are only visually amplified because the axes are truncated, sometimes showing only the top 2% of the domain. There are also some overly speculative statements, for example variant calling in microbes. The authors should edit the text to eliminate overstatements. 

<span style="color:red">
Not sure how this sentence can be deemed to be an overstaightment?
</span>

Second, Octopus is touted for somatic variant calling, but I am not convinced that what the authors show proves its worth here. They create and test with simulated data, something that continues to fall out of favor because of the massive corpus of actual cancer data that is now available. This testing evidently does not consider clonality, a biological factor that often dominates the analysis. There is likewise no obvious mention of handling sample impurity, a specimen-collection artifact that is hugely important for many types of cancer. The authors also test on tumor-only data, comparing Octopus to an unpublished caller. While there is some application for archival data, standardization of the tumor-normal paradigm makes this less relevant for newer data. 

<span style="color:red">
Explained in detail. Clonality/impurity is considered by virtue of testing a spectrum of allele frequencies.
</span>

The third criticism regards mathematical presentation of the Octopus algorithm, in which there seem to be mistakes and ambiguities, a few examples of which follow. On pp 11, the authors give p(h) from a coalescent model (ref 46) using what seem to be combinatorial coefficients (the "choose" notation) having 2nd index greater than the 1st, which is meaningless. Are these supposed to be fractions, but which are missing the solidus? Another is the overly casual way that symbol "g" (genotype) is used, e.g. to refer to disease (cancer), variant type (somatic, germline), sample (maternal, paternal). Indices (like 0 or 1) are sometimes wantonly added, as are absolute value bars to turn it into a number indicating ploidy of a sample. The bolded "g" indicates a set of genotypes, but sometimes so does plain "g". For example, the expression for "g union" on pp 12 uses set-builder notation that makes no sense ("h" is an element of "g" such that "g" is an element of "set g"). The
authors use other notation (like Iverson brackets), which many readers will find very difficult to follow and which can be expressed in a simpler form (e.g. using indices on the summation sign). Overall, I think the mathematical presentation needs attention.

<span style="color:red">
Okay
</span>