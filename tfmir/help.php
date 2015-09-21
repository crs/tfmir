<?php

?>
<!DOCTYPE html>
<html>

<head>
<title>TFmiR Help</title>
    <meta name="viewport" content="width=device-width, initial-scale=1.0">

    <link rel="stylesheet" href="faq.css">


</head>
<body>
<div class="pure-g">
<div class="pure-u-1-2"><img src="./img/overview.png" style="width:100%;"></div>
<div class="pure-u-1-2">	
<ul class="faq">
<li>
	<div class="q">What is this for?</div>
	<div class="a">
		<ol><li>For a list of up- and downregulated (1) mRNAs and (2) mRNAs from a differential analysis, 
		TFmiR builds a co-regulatory network based on parameters set (3).</li>
		<li>Database retrieval of predicted and experimentally validated interactions, based on confidence level</li>
		<li>Generate general and disease-specific co-regulatory networks</li>
		</ol>
	</div>
	</li>
	<li>
	<div class="q">What does the input look like?</div>
	<div class="a">
	The input consists of
		<ol> 
		<li> a list of mRNA identifiers with an indicator for their up- or downregulation (-1/1)</li>
		<li> (optional) a list of miRNA identifiers with an indicator for their up- or downregulation (-1/1)</li>
		<li> parameters that determine the p-value treshold for the overrepresentation analysis, the associated disease, and the confidence level (experimentally validated and/or predicted interactions)</li>
		</ol>
	</div>
	</li>
	<li>
	<div class="q">How exactly has a file to look like?</div>
	<div class="a">An input file (for mRNA, miRNA analogous) looks like:
	<pre>
AATK	-1
ABCB8	1
ABCG4	1
ABHD10	1
ABLIM1	-1
ABT1	1
etc...
	</pre>
	However, you can try the dataset we used for breast neoplasm:<br/>
	<a href="mirna.bc.txt" target="blank"><pre>mirna.bc.txt</pre></a>
	<a href="mrna.bc.txt" target="blank"><pre>mrna.bc.txt</pre></a><br />
	Alternatively, you can just click the &quot;Load example data&quot; at the start screen to load
	the dataset automatically.
	</div>
	</li>
	<li>
	<div class="q">Which databases are incorporated into TFmiR?</div>
	<div class="a">These databases are used:<br />
<!-- <img src="img/dbtable.png" alt="table of used databases"> -->
<table>
	<tr>
		<td>
			Interaction
		</td>
		<td>
			Databases
			(P/E) *
		</td>
		<td>
			Genes
		</td>
		<td>
			miRNAs
		</td>
		<td>
			Regulatory
			links
		</td>
		<td>
			Version
			/frozen date
		</td>
	</tr>
	<tr>
		<td>
			TF &#8594; gene
		</td>
		<td>
			<a title="Open database webpage" target="blank" href="http://www.gene-regulation.com/pub/databases.html">
			TRANSFAC</a>  (E) 
			<a title="Matys, V., Fricke, E., Geffers, R., Gößling, E., Haubrock, M., Hehl, R., Hornischer, K., Karas, D., Kel, A.E. and Kel-Margoulis, O.V. (2003) TRANSFAC®: transcriptional regulation, from patterns to profiles. Nucleic acids research, 31, 374-378" href="http://www.ncbi.nlm.nih.gov/pubmed?term=12520026" target="blank">[1]</a> 
		</td>
		<td>
			1279 
		</td>
		<td>
			--
		</td>
		<td>
			2943 
		</td>
		<td>
			V11.4
		</td>
	</tr>
	<tr>
	<td></td>
		<td>
			<a title="Open database webpage" target="blank" href="http://www.oreganno.org/">
			OregAnno
			</a> (E)
			<a title="Griffith, O.L., Montgomery, S.B., Bernier, B., Chu, B., Kasaian, K., Aerts, S., Mahony, S., Sleumer, M.C., Bilenky, M. and Haeussler, M. (2008) ORegAnno: an open-access community-driven resource for regulatory annotation. Nucleic acids research, 36, D107-D113" href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2239002/">[2]</a>
		</td>
		<td>
			1132
		</td>
		<td>
			--
		</td>
		<td>
			1083
		</td>
		<td>
			Nov 2010
		</td>
	</tr>
	<tr>
	<td></td>
		<td><a title="Open database webpage" target="blank" href="http://rulai.cshl.edu/TRED">
			TRED (P) 
			</a>
			<a title="Jiang, C., Xuan, Z., Zhao, F. and Zhang, M.Q. (2007) TRED: a transcriptional regulatory element database, new entries and other development. Nucleic acids research, 35, D137-D140" href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1899102/">[3]</a>
		</td>
		<td>
			3038
		</td>
		<td>
			--
		</td>
		<td>
			6462
		</td>
		<td>
			2007<a
			></a>
		</td>
	</tr>
	<tr>
		<td>
			TF &#8594; miRNA
		</td>
		<td>
		<a title="Open database webpage" target="blank" href="http://cmbi.bjmu.edu.cn/transmir">
			TransmiR
			</a> (E)
			<a title="Wang, J., Lu, M., Qiu, C. and Cui, Q. (2010) TransmiR: a transcription factor–microRNA regulation database. Nucleic acids research, 38, D119-D122" href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2808874/">[4]</a>
		</td>
		<td>
			158
		</td>
		<td>
			175
		</td>
		<td>
			567
		</td>
		<td>
			V1.2, Jan
			2013
		</td>
	</tr>
	<tr>
	<td></td>
		<td>
		<a title="Download dataset as published at BMC (.doc)" target="blank" href="http://www.biomedcentral.com/content/supplementary/1752-0509-4-90-s1.doc">
			PMID20584335
			</a> (E)
			<a title="Qiu, C., Wang, J., Yao, P., Wang, E. and Cui, Q. (2010) microRNA evolution in a human transcription factor and microRNA regulatory network. BMC systems biology, 4, 90" href="http://www.biomedcentral.com/1752-0509/4/90">[5]</a>
		</td>
		<td>
			58
		</td>
		<td>
			56
		</td>
		<td>
			102
		</td>
		<td>
			Apr 2009
		</td>
	</tr>
	<tr>
	<td></td>
		<td><a title="Open database webpage" target="blank" href="http://deepbase.sysu.edu.cn/chipbase/">
			ChipBase
			</a> (P)
			<a title="Yang, J.-H., Li, J.-H., Jiang, S., Zhou, H. and Qu, L.-H. (2013) ChIPBase: a database for decoding the transcriptional regulation of long non-coding RNA and microRNA genes from ChIP-Seq data. Nucleic acids research, 41, D177-D187" href="http://nar.oxfordjournals.org/content/41/D1/D177">[6]</a>
		</td>
		<td>
			119
		</td>
		<td>
			1380
		</td>
		<td>
			33087
		</td>
		<td>
			V1.1, Nov
			2012
		</td>
	</tr>
	<tr>
		<td>
			miRNA &#8594; gene
		</td>
		<td><a title="Open database webpage" target="blank" href="http://miRTarBase.mbc.nctu.edu.tw/">
			miRTarBase
			</a> (E)
			<a title="Hsu, S.-D., Lin, F.-M., Wu, W.-Y., Liang, C., Huang, W.-C., Chan, W.-L., Tsai, W.-T., Chen, G.-Z., Lee, C.-J. and Chiu, C.-M. (2010) miRTarBase: a database curates experimentally validated microRNA–target interactions. Nucleic acids research, gkq1107" href="http://www.ncbi.nlm.nih.gov/pubmed/21071411">[7]</a>
		</td>
		<td>
			2244
		</td>
		<td>
			551
		</td>
		<td>
			5640
		</td>
		<td>
			V4.5, Nov
			2013
		</td>
	</tr>
	<tr>
	<td></td>
		<td><a title="Open database webpage" target="blank" href="http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=tarbase/index">
			TarBase
			</a> (E)
			<a title="Sethupathy, P., Corda, B. and Hatzigeorgiou, A.G. (2006) TarBase: A comprehensive database of experimentally supported animal microRNA targets. Rna, 12, 192-197" href="http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1370898/">[8]</a>
		</td>
		<td>
			422
		</td>
		<td>
			79
		</td>
		<td>
			492
		</td>
		<td>
			V7.0
		</td>
	</tr>
	<tr><td></td>
		<td><a title="Open database webpage" target="blank" href="http://www.hsls.pitt.edu/obrc/index.php?page=URL1237998207">
			miRecords
			</a> (E)
			<a title="Xiao, F., Zuo, Z., Cai, G., Kang, S., Gao, X. and Li, T. (2009) miRecords: an integrated resource for microRNA–target interactions. Nucleic acids research, 37, D105-D110" href="http://nar.oxfordjournals.org/content/37/suppl_1/D105.abstract">[9]</a>
		</td>
		<td>
			543
		</td>
		<td>
			157
		</td>
		<td>
			780
		</td>
		<td>
			Mar 2009
		</td>
	</tr>
	<tr><td></td>
		<td><a title="Open database webpage" target="blank" href="http://starbase.sysu.edu.cn/">
			starBase
			</a> (P)
			<a title="Yang, J.-H., Li, J.-H., Shao, P., Zhou, H., Chen, Y.-Q. and Qu, L.-H. (2011) starBase: a database for exploring microRNA–mRNA interaction maps from Argonaute CLIP-Seq and Degradome-Seq data. Nucleic acids research, 39, D202-D209" href="http://www.ncbi.nlm.nih.gov/pubmed/21037263">[10]</a>
		</td>
		<td>
			5720
		</td>
		<td>
			249
		</td>
		<td>
			56051
		</td>
		<td>
			V2.0, Sept
			2013
		</td>
	</tr>
	<tr>
		<td>
			miRNA &#8594; miRNA
		</td>
		<td><a title="Open database webpage" target="blank" href="http://www.isical.ac.in/~bioinfo_miu/pmmr.php">
			PmmR</a>(P)
			<a title="Sengupta, D. and Bandyopadhyay, S. (2011) Participation of microRNAs in human interactome: extraction of microRNA–microRNA regulations. Molecular Biosystems, 7, 1966-1973." href="http://www.researchgate.net/publication/51042662_Participation_of_microRNAs_in_human_interactome_extraction_of_microRNA-microRNA_regulations">[11]</a>
			
		</td>
		<td>
			--
		</td>
		<td>
			312
		</td>
		<td>
			3846
		</td>
		<td>
			Mar 2011
		</td>
	</tr>
</table>
(P) means
predicted interactions and (E) means experimentally validated interactions.
	</div>
	</li>
	<li>
	<div class="q">I am still stuck!</div>
	<div class="a">Contact us for help:<br/>
	m.hamed@bioinformatik.uni-saarland.de<br/>
	christian.spaniol@bioinformatik.uni-saarland.de</div>
	</li>
</ul>
</div>
</div>

</body>

</html>