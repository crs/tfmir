<?php

require('parseFunctions.php');

session_start();

$interactions = array('tf-mirna', 'mirna-mirna', 'tf-gene', 'mirna-gene');

assert($_REQUEST['folder'] == session_id());

$folder = './uploads/' . $_REQUEST['folder'] . '/' . $_REQUEST['dataset'];
//echo $_REQUEST['folder'];
if (in_array($_REQUEST['dataset'], $interactions)) {
  echo getInteractionButtonSummary($folder);
} elseif ($_REQUEST['dataset'] == 'all' || $_REQUEST['dataset'] == 'disease') {
  echo getAllResultButtonSummaryDiv($folder);
} else {
  echo "You have no access here!";
}

function getInteractionButtonSummary($folder) {
 $summary = parseSummaryFile($folder . '/summary.txt');

 $template = 'Overlap significance (hyp.-geom. test): %e<br/>Simulation test p-val: %e';
 return sprintf($template, $summary['venn.pval.hypergeom'], $summary['venn.pval.simulation']);
}

function getAllResultButtonSummaryDiv($folder) {
  $summary = parseSummaryFile($folder . '/summary.txt');

  if ($_REQUEST['dataset'] == 'disease') {
  	$template = 'Graph with %d nodes and %d edges<br/>Network significance: %e';
  	return sprintf($template, $summary['graph.nodes.no'], $summary['graph.edges.no'], $summary['pval.hypergeom.disease.node']);	
  }
  else if ($_REQUEST['dataset'] == 'all') {
  	$template = 'Graph with %d nodes and %d edges';
  	return sprintf($template, $summary['graph.nodes.no'], $summary['graph.edges.no']);		
  }
}


?>