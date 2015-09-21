<?php

function validate_access() {
	if ((isset($_REQUEST['folder']) && isset($_REQUEST['disease'])) 
		&& preg_match("/[\w]{32}/", $_REQUEST['folder']) 
		&& (($_REQUEST['dataset'] == 'disease') || $_REQUEST['dataset'] == 'all')) {
		echo "ERROR! Corrupt session or dataset info";
		exit(1);
	}
}

function cmd_exec($cmd, &$stdout, &$stderr)
{
    $outfile = tempnam("/tmp/", "cmd");
    $errfile = tempnam("/tmp/", "cmd");
    $descriptorspec = array(
        0 => array("pipe", "r"),
        1 => array("file", $outfile, "w"),
        2 => array("file", $errfile, "w")
    );
    $proc = proc_open($cmd, $descriptorspec, $pipes);
    
    if (!is_resource($proc)) return 255;

    fclose($pipes[0]);    //Don't really want to give any input

    $exit = proc_close($proc);
    $stdout = file($outfile);
    $stderr = file($errfile);

//    unlink($outfile);
//    unlink($errfile);
    return $exit;
}

	function checkOption($select, $value) {
		if (isset($_COOKIE[$select]) && $_COOKIE[$select] == $value) {
			return 'selected="selected"';
		}
	}
	
	function checkRadio($select, $value) {
		if (isset($_COOKIE[$select]) && $_COOKIE[$select] == $value) {
			return 'checked="checked"';
		}
	
	}
	
	function checkBox($name) {
		if (isset($_COOKIE[$name]) && $_COOKIE[$name] == 'true') {
			return 'checked="checked"';
		}
	}
	
	function checkInput($name, $default) {
		if (isset($_COOKIE[$name])) {
			return $_COOKIE[$name];
		} else {
			return $default;
		}
	}

function folder_exists($file) {
	$complete_path = './uploads/' . session_id() . "/" . $file . "/";
	$files = glob($complete_path . "*");
	foreach ($files as $f) {
		if (is_file($f)) 
			return(basename($f));
			break;
	}
	return "None";
}

function deleteDir($dirPath) {
    if (! is_dir($dirPath)) {
        throw new InvalidArgumentException("$dirPath must be a directory");
    }
    if (substr($dirPath, strlen($dirPath) - 1, 1) != '/') {
        $dirPath .= '/';
    }
    $files = glob($dirPath . '*', GLOB_MARK);
    foreach ($files as $file) {
        if (is_dir($file)) {
            deleteDir($file);
        } else {
            unlink($file);
        }
    }
    rmdir($dirPath);
}	


?>
