<?php 
session_start();

include('functions.php');
//error_reporting(E_ALL);

$miRNA = $_COOKIE['miRNA'];
$mRNA = $_COOKIE['mRNA'];

$pval = $_POST['orapvalue'];
$disease = $_POST['disease'];
$output = $_POST['session'];
$evidence = $_POST['evidence'];


//if (!preg_match("/\w\.\w{1,2}/", $pval)) {
//	echo "pValue not valid";
//	$error = 1;
//}

//if (isset($error)) {
//	exit(1);
//}

$temp = getcwd() . '/uploads/' . $output;
mkdir($temp, 0777);

$mRNA = getcwd(). '/uploads/' . $output . '/mRNA/' . $mRNA;
$miRNA = getcwd(). '/uploads/' . $output . '/miRNA/' . $miRNA;

$output = getcwd() . '/uploads/' . $output;
$arguments = implode("','", array($mRNA, $miRNA, $pval, $evidence, $disease, $output));

$importPath = getcwd() . "/src/interface.R";

$diseasemotifs = $output . '/disease/motifs.txt';
$allmotifs = $output . '/all/motifs.txt';

if (file_exists($allmotifs)) {
	unlink($allmotifs);
} 

if (file_exists($diseasemotifs)) {
	unlink($diseasemotifs);
}

$in="source('" . $importPath . "'); start('".$arguments."');";

//echo $in;

$stderr = '';
$stdout = '';

echo "Submitted job, processing...";
echo "/usr/local/bin/R --vanilla -q -e \"$in\" >$temp/R3.Rout";
$exit_code = cmd_exec("/usr/local/bin/R --vanilla -q -e \"$in\" >$temp/R3.Rout", $stderr, $stdout);
echo "<br>Done, exit code: " . $exit_code;
//$exit_code = cmd_exec("ls -la", $stderr, $stdout);

//$exit_code = cmd_exec("R --vanilla -e \"print(2*2)\"", $a, $b); 
//print_r($a); print_r($b);
//echo "<pre>$in $exit_code</pre>";
//print_r($stderr);
//print_r($stdout);

?>