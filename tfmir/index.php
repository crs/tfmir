<?php
session_start();

error_reporting(0);

$line = date('Y-m-d H:i:s') . " - $_SERVER[REMOTE_ADDR] - $_SERVER[REMOTE_HOST] - $_SERVER[HTTP_REFERER]" ;
file_put_contents('./uploads/visitors.log', $line . PHP_EOL, FILE_APPEND);

$base = (!empty($_SERVER['HTTPS'])) ? "https://".$_SERVER['SERVER_NAME'].$_SERVER['REQUEST_URI'] : "http://".$_SERVER['SERVER_NAME'].$_SERVER['REQUEST_URI']; 
$base = "http://" . $_SERVER['HTTP_HOST'] . dirname($_SERVER['SCRIPT_NAME']);
$url = '//'.$_SERVER['HTTP_HOST'].dirname($_SERVER['PHP_SELF']);

$_SESSION['LAST_ACTIVITY'] = time();

?>

<!DOCTYPE html> 
<html>
<head>

<?php 
	include_once('./functions.php');
	//include_once('./cleanup.php'); 

?>

	<?php include('head.inc.html'); ?>
	  <!-- <link rel="stylesheet" href="./simplegrid.css" type="text/css"> -->

	<title>TFmiR - Data submission</title>
	<script type="text/javascript"> $(document).ready( checkResults );</script>
	<script type="text/javascript" src="./chosen_v1.2.0/chosen.jquery.min.js"></script>
	<link rel="stylesheet" href="./chosen_v1.2.0/chosen.min.css" type="text/css"> 

</head>
<body>

<div id="loading_dialog"></div>
<?php //echo realpath('.'); ?>
	<?php include('header.inc.html'); ?>

	<div id="ie" class="warning" style="display:none;">
	<h1>Sorry!</h1>
	<h2>Microsoft IE is not supported.</h2>
	<h4>It seems like you are using Internet Explorer, which is not supported yet.<br/>
	Please switch to <a href="https://www.google.de/chrome/browser/desktop/">Chrome</a>, <a href="http://download.mozilla.org">Firefox</a>, or <a href="http://www.apple.com/safari/">Safari</a>.</h4>
	</div>
	<script type="text/javascript">
	 if (msieversion()) {
	 	$('#ie').css('display','block');
	 }
	</script>
	<noscript>
	<div class="warning">
	<h1> Javascript required</h1>
<h4>This application needs Javascript enabled to function properly</h4>

</div>
	</noscript>
	<div id="pageframe" class="grid grid=pad">
		<div id="mainframe" class="col-9-12">
			<div class="panel" id="inputPanel">
				<div class="heading">Step 1: Input selection 
				<?php 
					if (!file_exists('./uploads/'.session_id())) {
						?><div class="exampleDataButton"><a href="#" onclick="loadExampleData(); return false;" class="starburst4" title="Load an example dataset for a functionality demonstration"><span><span><span><span><br />Load example data</span></span></span></span></a>
						</div><?php
					} else {
						?>
						<div class="exampleDataButton"><a href="#" onclick="clearSession(); return false;" class="starburst4" title="Load an example dataset for a functionality demonstration"><span><span><span><span><br />Clear session</span></span></span></span></a>
						</div><?php
					}?>
				</div>
				<table>
					<tr id="miRNArow" class="fileInputRow">
						<td class="desc">
							<label for="miRNAfileToUpload">miRNA</label>
						</td>
						<td class="filename" id="current_miRNA"><?php echo folder_exists('miRNA')?></td>
						<td>
							<form id="miRNA" enctype="multipart/form-data" method="post" action="<?php echo $url; ?>/upload.php?type=miRNA">
							<input type="hidden" name="type" value="miRNA"/><input id="mirnaDemo" type="hidden" name="demo" value="false" />
							<!-- <input id="miRNAfileToUpload-button" type="button" disabled="disabled" onclick="uploadFile('miRNA');" value="&lt;&lt;" title="Upload this file" />-->
							<input title="Input is a file where each line contains miRNA ID (e.g. 'hsa-mir-192a'), tabulator, and  1 or -1 for up- or downregulation" type="file" size="40" name="uploadedfile" id="miRNAfileToUpload" onchange="fileSelected(this); uploadFile('miRNA');"/></form>
						</td>						
					</tr>
					<tr id="mRNArow" class="fileInputRow">
						<td class="desc">
							<label for="mRNAfileToUpload">mRNA</label>
						</td>
						<td class="filename" id="current_mRNA"><?php echo folder_exists('mRNA')?></td>
						<td>
							<form id="mRNA" enctype="multipart/form-data" method="post" action="<?php echo $url; ?>/upload.php?type=mRNA">
							<input type="hidden" name="type" value="mRNA"/><input id="mrnaDemo" type="hidden" name="demo" value="false" />
							<!-- <input id="mRNAfileToUpload-button" type="button" disabled="disabled" onclick="uploadFile('mRNA');" value="&lt;&lt;" title="Upload this file"/>-->
							<input title="Input is a file where each line contains a gene symbol (e.g. 'GFAP', 'AATK'), tabulator, and  1 or -1 for up- or downregulation" or type="file" name="uploadedfile" id="mRNAfileToUpload" onchange="fileSelected(this); uploadFile('mRNA');"/></form>
						</td>
					</tr>
				</table>
			</div>
			<form method="post" id="execute" >
				<div class="panel">
					<div class="heading">Step 2: Configuration</div>
					<table>
						<tr>
							<td class="desc">p-Value treshold</td><td><input class="chosen_border" title="Enter a treshold between 0 and 1" name="orapvalue" type="text" size="6" value="<?php echo checkInput('orapvalue','0.05');?>" onchange="if (checkRange(0.0,1.0)) setCookie(this);"></td>
						</tr>
						<tr>
							<td class="desc">Related disease</td>
							<td>
								<select style="width:400px;" title="Choose a disease; use search field to narrow options" id="disease" name="disease" onchange="setCookie(this);">
								<option <?php echo checkOption('disease', ''); ?> value="">No disease</option><!--
									<option <?php checkOption('disease', 'bc'); ?> value="bc">Breast Cancer</option>
									<option <?php checkOption('disease', 'alz'); ?> value="alz">Alzheimer's Disease</option>
									<option <?php checkOption('disease', 'Melanoma'); ?> value="Melanoma">Melanoma</option> -->
									<?php 
									include_once('diseaseFunctions.php');
									echo getDiseaseOptions('../backend/disease.txt'); ?>
								</select>
							</td>
						</tr>
						<tr>
							<td class="desc">
								Evidence
							</td>
							<td>
								<select class="chosen_border" id="evidence" name="evidence" title="Choose evidence level; either predicted interactions, experimentally validated interactions, or both" onchange="setCookie(this);">
									<option <?php echo checkOption('evidence', 'Experimental') ?> value="Experimental">Experimental</option>
									<option <?php echo checkOption('evidence', 'Predicted') ?> value="Predicted">Predicted</option>
									<option <?php echo checkOption('evidence', 'Both') ?> value="Both">Both</option>
								</select>
							</td>
						</tr>
					</table>
				</div>
				
				

				<div class="panel" id="submitForm">
					<div class="heading">Step 3: Go!</div>
					<input type="hidden" id="mRNA-filename" name="mRNA-filename" value="NA" />
					<input type="hidden" id="miRNA-filename" name="miRNA-filename" value="NA" />
					<input type="hidden" id="session" name="session" value="<?php echo session_id(); ?>"/>
					<div id="startButtonContainer">
						<div id="processingButton" title="Please choose input files to enable processing" class="resultButton processButton center inactive"><!-- <a href="#" id="start" type="button" onclick="startProcessing(); return false;" value="Start processing!" /> -->
							<h4>Start processing</h4> <!-- < /a> --></div>
						</div>

					</div>
				</form>
				<div class="hiddenResultPanel">
				<div class="panel" id="result-chooser">
					<div class="heading">Step 4: Review result sets</div>
<br />
						<div class="heading subhead" title="Click on the network to get detailed information"> (a) Created networks</div>
						<div class="networkrow">

						<div id="leftDummy" class="resultButton inactive floating invisible">
							
						</div>


						<div id='disease' title="Click to explore the disease specific interaction network" class="diseaseButton resultButton inactive floating">
							<h3>Disease specific network</h3>
							<img class="interactionResultButton" src="img/all.png" alt="Disease specific interaction network" />
							<!-- <h4>Combined results</h4> -->
						</div>
						<div id='all' title="Click to explore the full interaction network" class="resultButton inactive floating">
							<h3>Full interaction network</h3>
							<img class="interactionResultButton" src="img/all.png" alt="Complete interaction network" />
							<!-- <h4>Combined results</h4> -->
						</div>
						</div>


					<div id="resultButtonsBar">	
					<div class="heading subhead">(b) Interaction types</div>
					<div class="grid grid-pad">
						<div class="col-1-2">
						<div id="tf-gene" title="Click to explore the transcription factor -&gt; gene interactions" class="resultButton inactive floating">
							<img class="interactionResultButton" src="img/tf-gene.png" alt="TF -> Gene interaction" />
							<!-- <h4>TF -> Gene</h4> -->
						</div>
						</div>
						<div class="col-1-2">
						<div id="tf-mirna" title="Click to explore the transcription factor -&gt; miRNA interactions" class="resultButton inactive floating">
							<img class="interactionResultButton" src="img/tf-mirna.png" alt="TF -> Gene interaction" />
							<!--<h4>TF -> miRNA</h4>-->
						</div>
						</div>
						<div class="col-1-2">
						<div id="mirna-mirna" title="Click to explore the miRNA -&gt; miRNA interactions" class="resultButton inactive floating">
							<img class="interactionResultButton" src="img/mirna-mirna.png" alt="miRNA -> miRNA interaction" />
							<!--<h4>miRNA -> miRNA</h4>-->
						</div>
						</div>

							
						<div id="mirna-gene" title="Click to explore the miRNA -&gt; gene interactions" class="resultButton inactive floating">
							<img class="interactionResultButton" src="img/mirna-gene.png" alt="miRNA -> Gene interaction" />
						</div>
						
						</div>
					</div>
				</div>
			</div>
			</div>
			<div class="col-3-12">
			<div class="log">
				<div class="heading">Log<div id="meter" class="meter"><span id="progressSpan"></span></div></div>
				<div id="currentJob">
					<div id="progressNumber"></div>
				</div>
				<div id="entries">
				</div>
			</div>
			</div>
		</div>

		<div id="footer">
			<?php echo date("Y"); ?> Mohamed Hamed, Christian Spaniol, Maryam Nazarieh, &amp; Volkhard Helms, <a href="http://gepard.bioinformatik.uni-saarland.de/" target="_blank">Chair of Computational Biology</a>
		</div>
	</body>
	<script type="text/javascript">
	$(function() {
		$('#disease').chosen();

		if (!msieversion()) {
			$('#pageframe').css('visibility','visible');
		}
	});
	</script>
	</html>
