<?php

session_start();

$baseDir = realpath("." );

$possibleDatasets = array('mirna-mirna', 'tf-mirna', 'tf-gene', 'mirna-gene', 'all');

//if (!in_array($_REQUEST['dataset'], $possibleDatasets) || $_REQUEST['folder'] != session_id()) {
//  echo "Access not allowed!";
//  exit(0);
//}

$subfolder = $baseDir . '/uploads/'. session_id() . '/'. $_REQUEST['dataset'];

//echo $subfolder;

include('parseFunctions.php'); 
include('dataset.php'); 

// preloading
//$all = getTabFor('all');
//$mirnamirna = getTabFor('mirna-mirna');

$tabs = getTabFor($_REQUEST['dataset']);

?><!DOCTYPE html>
<html>
<head>
<?php include('head.inc.html'); ?>
<title>tfmir result page</title>
</head>
<body>
<?php include('header.inc.html'); ?>
<div id="dialog"><code></code></div>
<div class="loading">Loading page<br /><img src="./img/loader.gif" alt="Please stand by..." /></div>
<div id="pageframe">
	<div id="tabs-nested-left" class="tabs">
		<ul class="ui-tabs-nav">
		<li class="ui-state-active ui-state-default ui-corner-top ui-tabs-anchor"><a class="ui-tabs-anchor" href="viewResult.php?dataset=disease">Disease<img class="sidebarimgLrg" src="./img/all.png" alt="Show disease-specific interaction network" /></a></li>
			<li class="ui-state-active ui-state-default ui-corner-top ui-tabs-anchor"><a class="ui-tabs-anchor" href="viewResult.php?dataset=all">Complete<img class="sidebarimgLrg" src="./img/all.png" alt="Show complete interaction network" /></a></li>
			<hr />
			<li class="ui-state-default ui-corner-top ui-tabs-anchor"><a class="ui-tabs-anchor" href="viewResult.php?dataset=mirna-mirna" id="mirna-mirna-tablink">miRNA-miRNA<img height="70"class="sidebarimg" src="./img/mirna-mirna.png" alt="show miRNA to miRNA interactions" /></a></li>
			<li class="ui-state-default ui-corner-top ui-tabs-anchor"><a class="ui-tabs-anchor" href="viewResult.php?dataset=tf-mirna">TF-miRNA<img height="70" class="sidebarimg" src="./img/tf-mirna.png" alt="show transcription factor to miRNA interactions" /></a></li>
			<li class="ui-state-default ui-corner-top ui-tabs-anchor"><a class="ui-tabs-anchor" href="viewResult.php?dataset=tf-gene">TF-Gene<img height="70" class="sidebarimg" src="./img/tf-gene.png" alt="show transcription factor to gene interactions" /></a></li>
			<li class="ui-state-default ui-corner-top ui-tabs-anchor"><a class="ui-tabs-anchor" href="viewResult.php?dataset=mirna-gene">miRNA-Gene<img height="70" class="sidebarimg" src="./img/mirna-gene.png" alt="show miRNA to gene interactions" /></a></li>
			<hr />
			<li class="ui-state-default ui-corner-top ui-tabs-anchor"><a class="ui-tabs-anchor" href="download.php?id=<?php echo session_id(); ?>" title="Click here to download a zip file with all generated results">Download results<img height="70" class="sidebarimg" src="./img/archive.png" alt="Click to download a zip archive containing all results" /></a></li>
			<!-- <li><a href="#tf-gene"><img class="sidebarimg" src="./img/mirna-gene.png" alt="show mirna to gene interactions" /></a></li> -->
		</ul>
		<div id="tab"><?php echo $tabs; ?></div>
		<!--
		<div id="complete"><?php echo $all; ?></div>
		<div id="mirna-mirna" onload="reloadNetwork('#mirna-mirna-network');"><?php //echo $mirnamirna; ?></div>
		<div id="tf-mirna">This is the tf-mirna tab</div>
		<div id="tf-gene">This is the tf-gene tab</div> -->
		<!-- <div id="mirna-gene">This is the forbidden fruit</div> -->
	</div>

</div>

</body>
<script src="cytoscape.js-panzoom.js"></script>
<script>

$(function() { 

//function bodyOnLoad() {
	console.log("Creating inner tabs bars");
	$(".tabs-min").tabs();
	console.log("Creating tab sidebar");
	//$("#tabs-nested-left").tabs();
	console.log("Show page");

	$('.loading').css('display', 'none');
	$('#pageframe').css('visibility','visible');
	//$('#all-network').cytoscape('get').container();
	//$('#all-network').cytoscape('get').resize();
	//$('#all-network').cytoscape('get').layout();
	//$('#all-network').cytoscape('get').center();
	var height = $('#all-network').height;
	var width = $('#all-network').width;
	
});

function reloadNetwork(network) {
	console.log('Reload Network ' + network);
	//$(network).css( { 'height': height, 'width': width });
	$(network).cytoscape('get').resize();
	$(network).cytoscape('get').layout();
}

//$("#mirna-mirna-tablink").click(function () { reloadNetwork('#mirna-mirna-network'); } ); 

$(document).ready(function() {
	console.log('page ready');

	//$.get('uploads/<?php echo session_id() . '/' . $_REQUEST['dataset']; ?>/motifs.txt', function(data) {
	$.get('checkMotif.php?dataset=<?php echo $_REQUEST['dataset'];?>', function(data) {
		console.log(data);
		if (data == '1') {
			motifSearch('<?php echo $_REQUEST['dataset']?>');
		}
	});
});
</script>
</html>

