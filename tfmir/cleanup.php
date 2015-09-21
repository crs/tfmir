<?php
	// cleanup folders older than a day

$dir = "./uploads/";

/*** cycle through all files in the directory ***/
foreach (glob($dir."*") as $file) {

/*** if file is one week (604800 seconds) old then delete it ***/
if (filemtime($file) < time() - 604800) {
    deleteDir($file);
    }
}


?>
