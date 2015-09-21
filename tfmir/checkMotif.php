<?php

session_start();

if (($_REQUEST['dataset'] != 'disease') && ($_REQUEST['dataset'] != 'all')) {
	echo "ERROR:\nWrong query";
}

$folder = 'uploads/';
if (file_exists($folder . session_id() . '/' . $_REQUEST['dataset'] . '/motifs.txt')) {
	echo 1;
} else {
	echo 0;
}

?>