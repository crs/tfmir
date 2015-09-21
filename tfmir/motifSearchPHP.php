<?php
$headers = array('Content-Type' => 'application/json');

session_start();
error_reporting(E_ALL);
//echo 'tst';

include('functions.php');

validate_access();

chdir('../backend');

$cmd_bin = "/usr/local/bin/R --vanilla -e \"%s\"";
//$expression = "source('TFMir.R'); PlotFunctionalSimilarity('".$gene_ids."','../tfmir/uploads/".$_REQUEST['folder']."/".$_REQUEST['motif'].".png')";


//$path_base = "../tfmir/uploads/".$_REQUEST['folder']."/".$_REQUEST['dataset'];
$path_base = "../tfmir/uploads/".session_id()."/".$_REQUEST['dataset'];

$expression = "source('TFMir.R'); ExportMotifs('".$path_base."/res.txt','".$path_base."','both')";

if (!file_exists($path_base."/motifs.txt")) {
	$stderr = '';
	$stdout = '';
	$exit_code = cmd_exec(sprintf($cmd_bin, $expression), $stderr, $stdout);
}
//echo $exit_code;

echo "<pre>";
print_r($stderr);
print_r($stdout);
echo "</pre>";


$motif_handle = fopen($path_base . "/motifs.txt", 'r');

$motif_template = '{ "type" : "%s", "tf" : "%s", "mirna" : "%s", "gene" : "%s" }';

$motifs = array();
if ($motif_handle) {
	$header = fgets($motif_handle);
    while (($line = fgets($motif_handle)) !== false) {
    	$fields = explode("\t", $line);
    	//print_r($fields);
    	array_push($motifs, sprintf($motif_template, $fields[1], $fields[2], $fields[3], $fields[4]));
    }
}

$pb = '{ "Motifs" : ['. implode(',', $motifs) . ']}';

echo $pb;
?>