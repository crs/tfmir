<?php
session_start();
include_once('./functions.php');
//header("Content-type: application/xhtml+xml");
//header("Access-Control-Allow-Origin: *")

if (($_REQUEST['type'] != 'miRNA') && ($_REQUEST['type'] != 'mRNA')) {
	//var_dump($_REQUEST);
	echo "Filetype was not given by form. This warning is not supposed to show.";
	exit();
}

$target_path = "./uploads/" . session_id() ."/";

$old_umask = umask(0);
if (!file_exists($target_path)) mkdir($target_path, 0777);

$target_path = "./uploads/" . session_id() ."/".$_REQUEST['type']."/";
if (file_exists($target_path)) {
	deleteDir($target_path);
}
mkdir($target_path, 0777);

//var_dump($_FILES);
$target_path = $target_path . basename( $_FILES['uploadedfile']['name']);

if(move_uploaded_file($_FILES['uploadedfile']['tmp_name'], $target_path)) {
		//setcookie($_REQUEST['type'], "../tf/uploads/" . session_id() ."/".$_REQUEST['type']."/".basename($_FILES['uploadedfile']['name']));
		setcookie($_REQUEST['type'], basename($_FILES['uploadedfile']['name']));
    echo " The file ".  basename( $_FILES['uploadedfile']['name']).
    " has been uploaded";
} else{
    echo "There was an error uploading the file ".$target_path.", please try again!";
}


?>
</html>
