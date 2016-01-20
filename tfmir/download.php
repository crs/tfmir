<?php

    require_once('Zip.php');

    error_reporting(4);

    if (isset($_GET['id'])) {
	    session_start($_GET['id']);
	} else {
	    session_start();
	}
    $filename = 'tfmir_' . substr(md5(session_id()),0,10);

    $source = './uploads/' . session_id();
    $target = './uploads/' . $filename .'.zip';
    
    if(file_exists($target)){
      unlink($target);
    }

    Zip($source, $target);


    header('Content-Type: application/zip');
    header("Content-Disposition: attachment; filename='". $filename . ".zip'");
    header('Content-Length: ' . filesize($filename));
    header("Location: " . $target);

?>
