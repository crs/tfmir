<?php

require('mappings.php');

function addGenecardLink($id) {
  return "<a href=\"http://www.genecards.org/cgi-bin/carddisp.pl?gene=".$id."\">".$id."</a>";
}

function urlCheck($text, $imgString = null) {
  // The Regular Expression filter
  //$reg_exUrl = "/(http|https|ftp|ftps)\:\/\/[a-zA-Z0-9\-\.]+\.[a-zA-Z]{2,3}(\/\S*)?/";
  $reg_exUrl = "/(http|https|ftp|ftps)\:\/\/[a-zA-Z0-9\-\.]+\.[a-zA-Z]{2,3}(\/(\S*\W*)*)?/";
  // Check if there is a url in the text
  if(preg_match($reg_exUrl, $text, $url)) {
    // make the urls hyper links
    if ($imgString == null) {
        return preg_replace($reg_exUrl, "<a href=\"{$url[0]}\">{$url[0]}</a>", $text);
    } else {
        $urlWithImage = "<a href=\"{$url[0]}\">".$imgString."</a>";
        return preg_replace($reg_exUrl, $urlWithImage, $text);
    }
    
  } else {
    // if no urls in the text just return the text
    return htmlspecialchars($text);
  }
}

//echo realpath('./uploads/test/mirna-mirna/mirna.ora.txt');
//echo parseTableFromCSV('./uploads/test/mirna-mirna/mirna.ora.txt', "oraTableStyle");

function align($str, $where)  { return '<div class="'.$where.'">'.$str.'</div>'; }
function lalign($str)         { return align($str, 'left'); }
function ralign($str)         { return align($str, 'right'); }
function calign($str)         { return align($str, 'center'); }

/**
 * @param $cell_no which cell from 1..n
 * @param $cell content of the cell
 * @param $type string of table header
 */
function oraTableStyle($cell_no, $cell, $type) {
  $pb = '';
  
  if (is_numeric($cell)) {
    $pb .= ralign(formatNumeric($cell));
  } elseif (oraCellTypes($cell_no) == 'array_string') {
      $pb .= lalign(lineBreak($cell));
  } else {
    $pb .= lalign($cell);
  } 
  
  
  return $pb;
}

function formatNumeric($number) {
  return number_format(round($number,2),2);
}

//function getResType($cell_no) {
//  $types = array('string', 'string', 'string', 'string', 'string', 'array_string', 'array_string', 'numeric', 'array_string', 'array_string', 'numeric', 'bool', 'bool');
//  return $types[$cell_no];
//}

function lineBreak($cell, $sep = ',') {
  $fragments = array_unique(explode($sep,$cell));
  return implode('<br />', $fragments);
}


function getDialogTableCell($cell, $id, $title, $sep = ",") {
  $pb = '';
  $diseases = array_unique(explode($sep,$cell));
  sort($diseases);
  $pb .= "<a href=\"#\" onclick=\"$('#".$id."').dialog({ title : '".$title."', position : {my : 'right bottom', at : 'right bottom', of: this} });return false;\">".count($diseases) . " hits</a>";

  $pb .= '<div id="'.$id.'" style="display:none;">';
  foreach ($diseases as $disease) { $pb .= $disease . '<br/>';}
  $pb .= '</div>';
  return $pb;
}


function defaultTableStyle($id, $cell, $title) {
  //print $cell;
  if (is_numeric($cell)) {
    if (strpos($cell, '.')) {
        return ralign(formatNumeric($cell));
    } else {
      return ralign($cell);
    }
  } elseif ((strpos($cell, ',') !== false) && strpos($title, 'Term') === false) {
    //$fragments = explode(',', $cell);
    //sort($fragments);
    //return lalign(implode('<br />', array_unique($fragments)));
    global $currentRegulator; global $currentTarget;
    //print $title;
    return getDialogTableCell($cell, $id, $title);
  } elseif ($cell == 'TRUE') {
    return calign('<img src="check.png" alt="association with disease" style="height:10pt;"/>');
  } elseif ($cell == 'FALSE'){
    return calign('<img src="x.png" alt="no association with disease" style="height:10pt;"/>');
  } else {
    return lalign($cell);
  }
  /*
  $type = getResType($cell_no);
  
  if ($type == 'string') {
    return(lalign($cell));
  } elseif($type == 'numeric') {
    return(ralign(formatNumeric($cell)));
  } elseif($type == 'array_string') {
    return(lalign(lineBreak($cell)));
  }
  return ($cell);*/
}


function getAligner($type) {
  if (strpos($type,'string') !== FALSE) {
    return 'lalign';
  } if (strpos($type, 'bool') !== FALSE) {
    return 'calign';
  } else {
    return 'ralign';
  }
}

function oraCellTypes($cell_no) {
  $types = array('string', 'string', 'int', 'array_string', 'numeric','numeric','numeric','numeric');
  return $types[$cell_no];
}

function oraCellTypeAligner($cell_no) {
  return getAligner(oraCellTypes($cell_no));
}

function getAlignmentClass($cell_no) {
  if (strpos(oraCellTypes($cell_no),'string') !== false) {
    return 'left';
  } else {
    return 'right';
  }
}

function resTableHead($head_no, $cell) {
  $pb = '<th>';
  
  /*
  switch ($head_no) {
    case 5: $pb .= 'Regulator Disease'; break;
    case 6: $pb .= 'Regulator Function'; break;
    case 7: $pb .= 'Regulator DSW'; break;
    case 8: $pb .= 'Target Function'; break;
    case 9: $pb .= 'Target Disease'; break;
    case 10: $pb .= 'Target DSW'; break;
    case 11: $pb .= 'Regulator in Disease'; break;
    case 12: $pb .= 'Target in Disease'; break;
    default: $pb .= ucwords($cell); break;
  }*/
  $pb .= getValue($cell);
  $pb .= '</th>';
  return $pb;
}



function defaultTableHead($head_no, $cell, $headKey) {
  $type = getKeyType($headKey);
  $aligner = getAligner($type);
  return "<th>" . $aligner(getValue($cell)) . "</th>";
}

function parseTableContentFromCSV($filename, $tableTypeFunction = 'defaultTableStyle', $tableTypeHead = 'defaultTableHead') {
  $script_src = '<script> $(document).ready( function() {';
  // final table string
  $pb = '';
  
  if (!file_exists($filename) || ($f = fopen($filename, "r")) === false) {
    return "Not available for this type of interaction";
  };
  
  // get table header
  $first = fgetcsv($f, 0, "\t");
  $pb .= '<thead><tr>';
  $head_no = 0;
  
  $header_cells = array();
  
  foreach ($first as $cell) {
    //$aligner = oraCellTypeAligner($head_no);
    //$pb .= '<th class="'.getAlignmentClass($head_no).'">' . $aligner(htmlspecialchars($cell)) . '</th>';
    
    array_push($header_cells, $cell);
    if (strpos($cell, 'genecard') === false) {
      $pb .= $tableTypeHead($head_no, $cell, $header_cells[$head_no]);
    } 
    $head_no += 1;
    
  }
  
  //print_r($header_cells);
  $pb .= '</tr></thead>';
  $pb .= '<tbody>';
  
  $line_no = 0;
  while (($line = fgetcsv($f, 0, "\t")) !== false) {
    $pb .= "<tr>";
    $cell_no = 0;
    foreach ($line as $cell) {
      
      $cell_id = substr(md5($filename),-4).$line_no.$cell_no;//$header_cells[$cell_no];

      if (preg_match('/^category$/', $header_cells[$cell_no])) {
        $interactions = explode('-', $cell);
        $regulator = $interactions[0];
        $target = $interactions[1];
      }
      //echo preg_match("/^regulator$/",$header_cells[$cell_no]);
      if ((preg_match("/^regulator$/",$header_cells[$cell_no])) && (!strstr($cell,'mir'))) {
        $currentRegulator = $cell;
        $cell = addGenecardLink($cell);
      }
      
      if ((preg_match("/^target$/", $header_cells[$cell_no])) && ((!strstr($cell,'mir')))) {
        $currentTarget = $cell;
        $cell = addGenecardLink($cell);
      }
      
      if (strpos($header_cells[$cell_no], 'genecard') !== false) {
        $cell_no += 1;
        continue;
      }
      //if (strpos($header_cells[$cell_no], 'genecard') !== false) {
      //  $cell = urlCheck($cell, '<img src="externalLink.png" />Genecard');
      //  $script_src .= "$('".md5($filename).$line_no."regulator').html(<a href=\"\"> + $('".md5($filename).$line_no."regulator').html() + </a>);\n";
      //}

      //$pb .= "<td id=\"".$cell_id."\">";      
      $pb .= "<td>";
      $pb .= $tableTypeFunction($cell_id, $cell, getValue($header_cells[$cell_no]));
      $pb .= "</td>";
      $cell_no += 1;
    }
    $pb .= "</tr>\n";
    $line_no += 1;
  }
  return($pb);// . "<tr style=\"display:none;\"><td>". $script_src." });</script><td></tr>");
}


function parseSummaryFile($file) {
  //$summary = parse_ini_file($file);
  require_once('INI-class.php');

  if (file_exists($file)) {
    $summary = new parseINI($file);
  }
  
  if ($summary->ini_array != null) {
    return array_pop($summary->ini_array);  
  } else {
    return false;
  }
  
}

?>
