<?php 

function endsWith($haystack, $needle) {
    // search forward starting from end minus needle length characters
    return $needle === "" || (($temp = strlen($haystack) - strlen($needle)) >= 0 && strpos($haystack, $needle, $temp) !== FALSE);
}

function checkOmitted($needle) {

    $omitted = array(
        "program.stderr",
        "program.stdout",
        "R3.Rout",
        "all.Rout",
        "finished.txt",
        "finishedMotifs.txt");

    foreach($omitted as $key) {
        if (endsWith($needle, $key)) {
            return true;
        }
    }
    return (false);
}


function Zip($source, $destination)
{
    $filenames = "";

    if (!extension_loaded('zip') || !file_exists($source)) {
        return false;
    }

    $zip = new ZipArchive();
    if (!$zip->open($destination, ZIPARCHIVE::CREATE)) {
        return false;
    }

    $source = str_replace('\\', '/', realpath($source));

    if (is_dir($source) === true)
    {
        $files = new RecursiveIteratorIterator(new RecursiveDirectoryIterator($source), RecursiveIteratorIterator::SELF_FIRST);

        foreach ($files as $file)
        {
            $file = str_replace('\\', '/', $file);

            if (checkOmitted($file)) {
                $filenames .= $file . " removed\n";
                continue;
            }

            // Ignore "." and ".." folders
            if( in_array(substr($file, strrpos($file, '/')+1), array('.', '..')) )
                continue;

            $file = realpath($file);

            if (is_dir($file) === true)
            {
                $zip->addEmptyDir(str_replace($source . '/', '', $file . '/'));
            }
            else if (is_file($file) === true)
            {
                $zip->addFromString(str_replace($source . '/', '', $file), file_get_contents($file));
            }
        }
    }
    else if (is_file($source) === true)
    {
        if (!checkOmitted($source))
            $zip->addFromString(basename($source), file_get_contents($source));
    }

    //$zip->addFromString('ziplog.txt', $filenames);
    return $zip->close();
}
?>
