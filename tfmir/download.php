<?php

    require_once('Zip.php');

    error_reporting(4);
    session_start();

    $filename = 'tfmir_' . substr(md5(session_id()),0,5);

    $source = './uploads/' . session_id();
    $target = './uploads/' . $filename .'.zip';
    
    Zip($source, $target);


    header('Content-Type: application/zip');
    header("Content-Disposition: attachment; filename='". $filename . ".zip'");
    header('Content-Length: ' . filesize($filename));
    header("Location: " . $target);

?>