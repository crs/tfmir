<?php

function urlCheck($text) {
  // The Regular Expression filter
$reg_exUrl = "/(http|https|ftp|ftps)\:\/\/[a-zA-Z0-9\-\.]+\.[a-zA-Z]{2,3}(\/\S*)?/";

// Check if there is a url in the text
if(preg_match($reg_exUrl, $text, $url)) {

       // make the urls hyper links
       return preg_replace($reg_exUrl, "<a href=\"{$url[0]}\">{$url[0]}</a>", $text);

} else {

       // if no urls in the text just return the text
       return htmlspecialchars($text);
}
}


function getmiRNATable($filename) {
  return '<table id="miRNATable">'.getResultTableLines($filename).'</table>';
}

function getDialogTableCell($cell, $id, $title, $sep = ",") {
  $pb = '';
  $diseases = explode($sep,$cell);
  $pb .= "<a href=\"#\" onclick=\"$('#".$id."').dialog({ title : '".$title."', position : {my : 'right bottom', at : 'right bottom', of: this} });return false;\">".count($diseases) . "</a>";

  $pb .= '<div id="'.$id.'" style="display:none;">';
  foreach ($diseases as $disease) { $pb .= $disease . '<br/>';}
  $pb .= '</div>';
  return $pb;
}

function getResultTableLines($filename) {
$pb = '';
$f = fopen($filename, "r");
$first = fgetcsv($f, 0, "\t");
$pb .= '<thead>';
foreach ($first as $cell) {
  $pb .= '<th>' . htmlspecialchars($cell) . '</th>';
}
$pb .= '</thead><tbody>';

$line_no = 0;
while (($line = fgetcsv($f, 0, "\t")) !== false) {
        $pb .= "<tr>";
        $cell_no = 0;
        foreach ($line as $cell) {
                if ($cell_no == 0) {
                  $current_gene = $cell;
                }
                if ($cell_no == 1) {
                  $current_mirna = $cell;
                }
                $pb .= "<td>";
                if ($cell_no == 6) {
                  $pb .= getDialogTableCell($cell, 'miRNA-disease'.$line_no.$cell_no, $current_mirna. " diseases", ';');
                } elseif ($cell_no == 12) {
                  $pb .= getDialogTableCell($cell, 'gene-disease'.$line_no.$cell_no, $current_gene. " diseases");
                } elseif ($cell_no == 16) {
                    $pb .= getDialogTableCell($cell, 'mir-association-disease'.$line_no.$cell_no, $current_mirna. " associated diseases");
                } elseif ($cell_no == 13) {
                    $pb .= '<a href="'.$cell.'" title="Open item in David database">Open David</a>';                  
                } elseif ($cell_no == 14) {
                    $pb .= '<a href="'.$cell.'" title="Open item in Genecard database">Open GeneCard</a>';
                } else {
                  $pb .= urlCheck($cell);
                }
                $pb .= "</td>";
                $cell_no += 1;
        }
        $pb .= "</tr>\n";
}
fclose($f);
return $pb . '</tbody>';
}

echo '<!-- ' . $_REQUEST['miRNATableFile'] . '-->';
if (file_exists($_REQUEST['miRNATableFile'])) {
  //echo 'The file is available.';
}
echo getmiRNATable($_REQUEST['miRNATableFile']);

?>
