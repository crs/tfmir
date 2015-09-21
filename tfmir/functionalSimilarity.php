<?php 

require('./functions.php');
error_reporting(E_ALL);
//echo shell_exec("/Library/Frameworks/R.framework/Resources/bin/R --vanilla -e 'print(version)'");




$stdout = "";
$stderr = "";

//cmd_exec("export PATH=/usr/lib/:/usr/local/lib & /Library/Frameworks/R.framework/Resources/bin/R --vanilla -e 'print(version)'", $stdout, $stderr);
//cmd_exec("export PATH=/Users/chris/Dropbox/Scripts & lottery.py", $stdout, $stderr);
//cmd_exec("sh --login -c 'R --vanilla -e \"print(version)\"'", $stdout, $stderr);
$saved = getenv("PATH") ."\n";

//cmd_exec("PATH=/Applications/XAMPP/xamppfiles/lib/:/usr/local/lib:/usr/bin/:/Library/Frameworks/R.framework/Versions/3.1/Resources/lib/ LD_LIBRARY_PATH=/Applications/XAMPP/xamppfiles/lib/:/usr/lib/:/usr/local/lib:/Library/Frameworks/R.framework/Versions/3.1/Resources/lib/ R --vanilla -e 'print(version)'", $stdout, $stderr);

$map = file_get_contents('./uploads/' . $_REQUEST['folder'] . '/genemap.txt');
//cmd_exec("DYLD_LIBRARY_PATH=/usr/lib:\$DYLD_LIBRARY_PATH R --vanilla --verbose -e 'print(tempdir())'", $stdout, $stderr);
//cmd_exec("DYLD_LIBRARY_PATH=/Library/Frameworks/R.framework/Versions/3.1/Resources/library:/usr/lib:\$DYLD_LIBRARY_PATH R -f '../../backend/TFMir.R'", $stdout, $stderr);
//cmd_exec("ls -la", $stdout, $stderr);
//print_r($stderr);
//print_r($stdout);

//echo $map;

$array = explode("\n", $map);


$new_array = array();


foreach ($array as $tuple) {
    $nums = explode("\t",$tuple);
    if (isset($nums[0]) && isset($nums[1])) {
        $new_array[$nums[0]] = $nums[1];
    }
}

$genes = explode(',', $_REQUEST['genes']);

$gene_ids = array();
foreach($genes as $gene_name) {
    if (array_key_exists($gene_name, $new_array)) {
         //$new_array[$gene_name] . "\n";    
         array_push($gene_ids, $new_array[$gene_name]);
    }
}
print_r($gene_ids);
$gene_ids = implode(',', $gene_ids);

//cmd_exec("DYLD_LIBRARY_PATH=/usr/lib:\$DYLD_LIBRARY_PATH R --vanilla --verbose -e 'print(tempdir())'", $stderr, $stdout);
chdir('../backend');

$cmd_bin = "/usr/local/bin/R --vanilla -e \"%s\"";
//$expression = "source('TFMir.R'); PlotFunctionalSimilarity('".$gene_ids."','../tfmir/uploads/".$_REQUEST['folder']."/".$_REQUEST['motif'].".png')";

$genes = implode(',', $genes);
$expression = "source('TFMir.R'); PlotFunctionalSimilarity('".$genes."','../tfmir/uploads/".$_REQUEST['folder']."/".$_REQUEST['type'].'-'.$_REQUEST['motif'].".png')";
$exit_code = cmd_exec(sprintf($cmd_bin, $expression), $stderr, $stdout);


//print_r($_REQUEST);
//print_r($genes);   
//print_r($array);
//print_r($new_array);

print($exit_code ."\n");
print_r($stderr);
print_r($stdout);
?>