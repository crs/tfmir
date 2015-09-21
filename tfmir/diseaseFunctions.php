<?php


function getDiseaseOptions($filename) {
  
  $f = fopen($filename, 'r');
  $pb = '';
  if ($f) {
    while (($line = fgets($f)) !== false) {
    	$line = trim($line);

    	if (mb_detect_encoding($line, 'UTF-8', true)) {
        	$pb .= '<option '. checkOption('disease', $line).' value="'.$line.'">' . $line . '</option>';
    	}
    }
  }
  return $pb;
}

?>