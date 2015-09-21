<!DOCTYPE html>
<html>
<head>
<?php include("head.inc.html");?>
<title>
<?php echo $_GET['tf']; ?>-<?php echo $_GET['mirna']; ?>-<?php echo $_GET['gene']; ?> Motif
</title>
</head>
<body>

<script type="text/javascript">;

var tf = '<?php echo $_GET['tf'];?>';
var mirna = '<?php echo $_GET['mirna'];?>';
var gene = '<?php echo $_GET['gene'];?>';

var motif = tf+'-'+mirna+'-'+gene;
var co = '<?php echo $_REQUEST['type']; ?>';
</script>
<h3>
<?php

	if ($_REQUEST['type'] == 'cotargeted') {
		$type = 'Co-targeted';
		$co = $_REQUEST['type'];
	} else {
		$type = 'Co-regulated';
		$co = $_REQUEST['type'];
	}
	echo $type;
?> subnetwork for TF: <?php echo $_GET['tf']; ?>, miRNA: <?php echo $_GET['mirna']; ?>, Gene: <?php echo $_GET['gene']; ?></h3>

<div><a href="#" onclick="exportImage(); return false;" title="Export current view as PNG image">Export network as image</a></div>
<div id="motifNetwork" class="cytoPanel" style="height:600px;width:1000px;"></div>


<div id="funcSimImage">
<?php 
	$motif = $_REQUEST['tf']."-".$_REQUEST['mirna']."-".$_REQUEST['gene'];

	$image_filename = "./uploads/".$_REQUEST['folder']."/".$co.'-'.$motif.".png";
	//echo $image_filename;
	if (file_exists($image_filename)) {
		//echo '<div><a href="#" onclick="calculateFunctionalSimilarity(\'' . $_REQUEST['folder'] .'\',pb,motif,co);return false;">Re-calculate functional similarity</a><span id="similarityLoader"><img height="15" src="img/loading.gif">calculating... please do not reload this page</span></div>';
		echo "<img src=\"".$image_filename."\" />";
	} else {
		echo '<div><a id="calculationButtonFS" href="#" title="Attention: this will take some time" onclick="calculateFunctionalSimilarity(\'' . $_REQUEST['folder'] .'\',pb,motif,co);return false;">Calculate functional similarity</a><span id="similarityLoader"><img height="15" src="img/loading.gif">calculating... please do not reload this page</span></div>';
	}	
?>
</div>
<!-- <div><a href="#" onclick="calculateFunctionalSimilarity('<?php echo $_REQUEST['folder'] ?>',pb,motif);return false;">Calculate functional similarity</a></div> -->
<!-- <script src="jquery-1.11.1.js" type="text/javascript"></script>
<script src="build/cytoscape.min.js" type="text/javascript"></script>-->
<script type="text/javascript">
var elesJson = <?php require('createCytoscapeNetwork.php'); echo createJSON('./uploads/' . $_REQUEST['folder'] . '/'.$_REQUEST['dataset'].'/res.txt') .";" ?>
</script>
<script src="motifNetwork.js" type="text/javascript">
</script>
<script src="functionalSimilarity.js.php" type="text/javascript"></script>
<script type="text/javascript">

	//$(document).ready( function () {

    var cy = $('#motifNetwork').cytoscape('get');

	var selector = 'node#'+tf+',node#'+mirna+',node#'+gene;
    
    var elements = cy.elements('node#'+tf+',node#'+mirna+',node#'+gene);
    //elements.css('-webkit-box-shadow','0px 0px 15px 5px rgba(255, 255, 190, .75)');
	//elements.css('-moz-box-shadow', '0px 0px 15px 5px rgba(255, 255, 190, .75)');
	//elements.css('box-shadow','0px 0px 15px 5px rgba(255, 255, 190, .75)');
	elements.css('border-width','5px');
    //elements = elements.add(elements.neighborhood());
    //elements = elements.add(cy.nodes("[source='"+tf+"'],[source='"+mirna+"']"]));
    //elements = elements.add(cy.edges("[source='"+tf+"'],[source='"+mirna+"']"));

    var type = '<?php echo $type; ?>'

    
    if (type == 'Co-targeted') {
    	//elements = elements.not(cy.elements('node#'+tf).outgoers().intersect(cy.elements('node#'+mirna).outgoers())).remove();
    	//elements = elements.add(cy.elements("[source='"+tf+"']").intersect(cy.elements("[source='"+mirna+"']")));  	
 		//elements = elements.add(cy.elements("[source='"+tf+"']").targets().intersect(cy.elements("[source='"+mirna+"']").targets()));  	
 		//elements = elements.add(cy.elements("edge[source='"+tf+"'],node[source='"+tf+"']"));
 		tf_outgoers = cy.elements('node#'+tf).outgoers();
 		mirna_outgoers =  cy.elements('node#'+mirna).outgoers();
 		common = tf_outgoers.intersect(mirna_outgoers);	
 		elements = elements.add(common);
 		elements = elements.add(common.incomers('node#'+tf+',node#'+mirna+',edge[source="'+tf+'"],edge[source="'+mirna+'"]'));
    } else {
    	elements = elements.add(cy.elements('node#'+tf+',node#'+mirna).outgoers());
    }
    //cy.elements().not(elements).select();
    cy.elements().not(elements).remove();

    var pb = '';
    cy.elements().forEach( function (ele, i, eles) {
    	if (ele.isNode() && ele.data('type') != 'mirna') {
    		pb = pb + ele.id() + ','; //console.log(ele.id());
    	}
    });

    //$(function() { $(document).tooltip(); });
</script>
<script type="text/javascript" src="cytoscape.js-panzoom.js"></script>
<script type="text/javascript">$(document).ready(function () { $('.cytoPanel').cytoscape('get').panzoom(); });</script>
</body>

</html>