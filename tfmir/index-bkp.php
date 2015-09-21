<?php
	session_start();

//	if (isset($_SESSION['LAST_ACTIVITY']) && (time() - $_SESSION['LAST_ACTIVITY'] > 1800)) {
//    // last request was more than 30 minutes ago
//    session_unset();     // unset $_SESSION variable for the run-time
//    session_destroy();   // destroy session data in storage
//}
$_SESSION['LAST_ACTIVITY'] = time();

	include_once('./functions.php');
	include_once('./cleanup.php');

	/*

	echo '<pre>';
	var_dump(session_id());
	var_dump($_COOKIE);
	folder_exists('miRNA');
	echo "</pre>";
	*/

?>
<html>
<head>
	<?php include('head.inc.html'); ?>
	<title>Transcription Factor Prediction Server</title>
	<script type="text/javascript"> $(document).ready( checkResults );</script>
</head>
<body >
<!--
		<div id="fileName"></div>
    <div id="fileSize"></div>
    <div id="fileType"></div>
		<div id="progressNumber"></div>
		-->

<!--
		<div id="fileName"></div>
    <div id="fileSize"></div>
    <div id="fileType"></div>
		<div id="progressNumber"></div>
		-->
<div id="header">
<div class="logoRight"><a href="http://www.zbi.uni-saarland.de/"><img src="img/logo-cbi.png" style="height:75%;"/></a></div><H2> Transcription Factor Prediction</H2>
</div>

<div class="log">
<div class="heading">Log<div id="meter" class="meter"><span id="progressSpan"></span></div></div>
	<div id="currentJob">
		<div id="progressNumber"></div>
	</div>
	<div id="entries">
	</div>
</div>

<div id="pageframe">
<div id="mainframe">
<div class="panel" id="inputPanel">
<div class="heading">Step 1: Input selection</div>
	<table>
	<tr id="miRNArow">

		<form id="miRNA" enctype="multipart/form-data" method="post" action="http://localhost/tf/upload.php?type=miRNA">
		<input type="hidden" name="type" value="miRNA"/>
		<td class="desc">
			<label for="miRNAfileToUpload">miRNA</label>
		</td>
		<td class="filename" id="current_miRNA"><?php echo folder_exists('miRNA')?></td>
		<td>
      <input id="miRNAbutton" type="button" disabled="disabled" onclick="uploadFile('miRNA');" value="&lt;&lt;" title="Upload this file" /><input type="file" name="uploadedfile" id="miRNAfileToUpload" onchange="fileSelected('miRNAbutton');"/>
    </td>
		</form>

	</tr>
	<tr id="mRNArow">

		<form id="mRNA" enctype="multipart/form-data" method="post" action="http://localhost/tf/upload.php?type=mRNA">
		<input type="hidden" name="type" value="mRNA"/>
			<td class="desc">
				<label for="mRNAfileToUpload">mRNA</label>
			</td>
	<td class="filename" id="current_mRNA"><?php echo folder_exists('mRNA')?></td>
	<td>
      <input id="mRNAbutton" type="button" disabled="disabled" onclick="uploadFile('mRNA');" value="&lt;&lt;" title="Upload this file"/><input type="file" name="uploadedfile" id="mRNAfileToUpload" onchange="fileSelected('mRNAbutton');"/>
    </td>
		</form>
		<form method="post" id="execute" >
		</tr>

</table>
</div>

<div class="panel">
		<div class="heading">Step 2: Configuration</div>
			<table>
				<tr>
				<td class="desc">p-Value treshold</td><td><input  title="Enter a treshold between 0 and 1" name="orapvalue" type="text" size="6" value="<?php checkInput('orapvalue','0.05');?>" onchange="if (checkRange(0.0,1.0)) setCookie();"></td>
			</tr>
	<tr>
		<td class="desc">Related disease</td>
		<td>
			<select style="width:400px;" id="disease" name="disease" onchange="setCookie();">
				<option <?php checkOption('disease', 'None'); ?> value="">No disease</option><!--
				<option <?php checkOption('disease', 'bc'); ?> value="bc">Breast Cancer</option>
				<option <?php checkOption('disease', 'alz'); ?> value="alz">Alzheimer's Disease</option>
				<option <?php checkOption('disease', 'Melanoma'); ?> value="Melanoma">Melanoma</option> -->
				<?php 
					include('diseaseFunctions.php');
					
					echo getDiseaseOptions('../backend/disease.txt');
				?>
			</select>
		</td>
	</tr>
	<tr>
			<td class="desc">
				Evidence
			</td>
			<td>
				<select name="evidence" onchange="setCookie();">
					<option <?php checkOption('evidence', 'Experimental') ?> value="Experimental">Experimental</option>
					<option <?php checkOption('evidence', 'Predicted') ?> value="Predicted">Predicted</option>
					<option <?php checkOption('evidence', 'Both') ?> value="Both">Both</option>
				</select>
			</td>
	</tr>
	</table>
</div>


	<!-- <input type="checkbox" name="testCheck" value="TRUE">   this box checked?-->
	<!--
	 <form id="form1" enctype="multipart/form-data" method="post" action="http://localhost/tf/upload.php">
    <div class="row">
      <label for="fileToUpload">Select a File to Upload</label><br />
      <input type="file" name="uploadedfile" id="fileToUpload" onchange="fileSelected();"/>
    </div>
    <div id="fileName"></div>
    <div id="fileSize"></div>
    <div id="fileType"></div>
    <div class="row">
      <input type="button" onclick="uploadFile();" value="Upload" />
    </div>

  </form>
	-->





 <!-- action="./cgi-bin/execute.cgi"> -->
<!--
<div class="panel">
	<div class="heading">ORA</div>
		<table>
			<tr>
				<td class="desc">Perform ORA for</td><td><input <?php checkBox('oragenes'); ?> name="oragenes" onchange="setCookie();" type="checkbox">Genes</td><td><input <?php checkBox('oramirna'); ?> onchange="setCookie();" name="oramirna" type="checkbox">miRNA</td>
			</tr>
			<tr>
				<td class="desc">p-Value treshold</tdclass="desc"><td><input  title="Enter a treshold between 0 and 1" name="orapvalue" type="text" size="6" value="<?php checkInput('orapvalue','0.05');?>" onchange="if (checkRange(0.0,1.0)) setCookie();"></td>
			</tr>
			<tr>
				<td class="desc">Multiple test correction</td>
					<td><input type="radio" <?php checkRadio('mtc','false') ?> onchange="setCookie();" name="mtc" value="false">None</td>
					<td><input type="radio"  <?php checkRadio('mtc','bh') ?> onchange="setCookie();" name="mtc" value="bh">BH</td>
				</tr>
				<tr>
					<td></td>
					<td><input type="radio"  <?php checkRadio('mtc','bonferroni') ?> onchange="setCookie();" name="mtc" value="bonferroni">Bonferroni</td>
					<td><input type="radio"  <?php checkRadio('mtc','fdr') ?> onchange="setCookie();" name="mtc" value="fdr">FDR</td>
			</tr>
			<tr>
				<td class="desc">Enrichment for mRNA/TF</td>
				<td colspan=2>
					<select multiple="multiple" size="3" title="Hold shift and click to select multiple data sources" name='mrnatf_enrichment'>
						<option value="go">GO</option>
						<option value="kegg">KEGG</option>
						<option value="hmdd">HMDD</option>
					</select>
				</td>

			</tr>
			<tr>
				<td class="desc">Enrichment for miRNA</td>
				<td colspan='2'>
					<select multiple="multiple" size="3" name='mirna_enrichment'>
						<option value="go">GO</option>
						<option value="kegg">KEGG</option>
						<option value="hmdd">HMDD</option>
					</select>
			</td>

			</tr>
	</table>
</div>
-->

<!--
<div class="panel">
	<div class="heading">Configuration</div>
	<table>
		<tr>
			<td class="desc">Related disease</td>
			<td colspan="2">
				<select name="disease" onchange="setCookie();">
					<option <?php checkOption('disease', 'None'); ?> value="None">No disease</option>
					<option <?php checkOption('disease', 'bc'); ?> value="bc">Breast Cancer</option>
					<option <?php checkOption('disease', 'alz'); ?> value="alz">Alzheimer's Disease</option>
				</select>
			</td>
		</tr>
		<tr>
			<td class="desc">Interaction Types</td>
			<td><input <?php checkBox('genes-tf'); ?> onchange="setCookie();" type="checkbox" name="genes-tf">TF-Genes</td>
			<td><input type="checkbox" name="mirna-mrna" onchange="setCookie();" <?php checkBox('mirna-mrna'); ?>/>miRNA-mRNA</td>
			<td><input onchange="setCookie();" <?php checkBox('tf-mirna'); ?> type="checkbox" name="tf-mirna">TF-miRNA</td>
		</tr>
		<tr>
			<td class="desc">Dependency test p-Value treshold</td>
			<td colspan="2"><input title="Enter a treshold between 0 and 1" onchange="if (checkRange(0.0,1.0)) setCookie();" type="text" name="deptestpval" size="6" value="<?php checkInput('deptestpval','0.05');?>"></td>
		</tr>
		<tr>
			<td class="desc">Confidence level</td>
			<td><input type="checkbox" onchange="setCookie();" <?php checkBox('expvalidated'); ?> name="expvalidated">Experimentally validated</td>
			<td><input onchange="setCookie();" <?php checkBox('predicted'); ?> type="checkbox" name="predicted">Predicted</td>
		</tr>
	</table>
</div >
-->
<div class="panel" id="submitForm">
	<div class="heading">Step 3: Go!</div>
	<input type="hidden" id="mRNA-filename" name="mRNA-filename" value="NA" />
	<input type="hidden" id="miRNA-filename" name="miRNA-filename" value="NA" />
	<input type="hidden" id="session" name="session" value="<?php echo session_id(); ?>"/>
	<div id="startButtonContainer">
		<div id="processingButton" class="resultButton processButton center" onClick="startProcessing(); return false;"><!-- <a href="#" id="start" type="button" onclick="startProcessing(); return false;" value="Start processing!" /> -->
		<h4>Start processing</h4> <!-- < /a> --></div>
	</div>
</form>
</div>
<div class="panel" id="result-chooser">
	<div class="heading">Step 4: Review result sets</div>
	
	
	
	<div id="resultButtonsBar">	
		<div id="tf-gene" class="resultButton inactive floating">
			<img class="interactionResultButton" src="img/tf-gene.png" alt="TF -> Gene interaction" />
			<!-- <h4>TF -> Gene</h4> -->
		</div>
		<div id="tf-mirna" class="resultButton inactive floating">
			<img class="interactionResultButton" src="img/tf-mirna.png" alt="TF -> Gene interaction" />
			<!--<h4>TF -> miRNA</h4>-->
		</div>
		<div id="mirna-mirna" class="resultButton inactive floating">
			<img class="interactionResultButton" src="img/mirna-mirna.png" alt="TF -> Gene interaction" />
		<!--<h4>miRNA -> miRNA</h4>-->
		</div>
	
	<!-- 
	===
	This is hidden because of reasons 
	===
		<div id="mirna-gene" class="resultButton inactive floating">
			<img class="interactionResultButton" src="img/mirna-gene.png" alt="TF -> Gene interaction" />
		</div>
	-->
		<div id='all' class="resultButton inactive floating">
			<img class="interactionResultButton" src="img/all.png" alt="TF -> Gene interaction" />
			<!-- <h4>Combined results</h4> -->
		</div>
	</div>
</div>
<div class="panel" id="results" style="display:none;">
	<div class="heading">Step 4: Result assessment</div>
	<div id='resultBox'></div>
</div>
</div>
</div>
<br />
<div id="footer">
 <?php echo date("Y"); ?> Mohamed Hamed &amp Christian Spaniol, <a href="http://gepard.bioinformatik.uni-saarland.de/" target="_blank">Chair of Computational Biology</a>
</div>
</body>
</html>
