<?php
$possibleDatasets = array('mirna-mirna', 'tf-mirna', 'tf-gene', 'mirna-gene', 'all');

//if (!in_array($_REQUEST['dataset'], $possibleDatasets) || $_REQUEST['folder'] != session_id()) {
//  echo "Access not allowed!";
//  exit(0);
//}

//$subfolder = $baseDir . '/uploads/'. session_id() . '/'. $_REQUEST['dataset'];

//echo $subfolder;


function getTabFor($dataset) {
  $baseDir = realpath("." );
  $subfolder = $baseDir . '/uploads/'. session_id() . '/';//. $_REQUEST['dataset'];
  
  $subfolder .= $dataset;
  $network = $dataset . '-network';
  $summary = $dataset . '-summary';
  $interactions = $dataset . '-interactions';
  $ora = $dataset . '-ora';
  $oraTable = $dataset . '-oraTable';
  $resTable = $dataset . '-resTable';

  $motif = $dataset . '-motif';
  $motifTable = $dataset . '-motifTable';

  if (!file_exists($subfolder . '/summary.txt')) {
    return 'This data is not available for the given dataset.';
  }

  $pb = "<script> $(function() { /*$(\"#tabs-min\").tabs()*/; $('#$motifTable').DataTable(); $('#".$oraTable."').DataTable(); $('#geneORATable').DataTable(); $('#".$resTable."').DataTable(); }); </script>";
  $pb .= "<script src=\"cytoscape.php?folder=". session_id() . "&dataset=".$dataset."\" type=\"text/javascript\"> </script>";

  $pb .= "<div class=\"tabs-min\"><ul><li><a href=\"#".$network."\" onclick=\"console.log('Invoke resize');$('#".$network."').resize(); console.log('Resize complete');\">Network</a>";

  //$pb .= "</li><li><a href=\"#".$summary."\">Summary</a>";

  $pb .= "</li><li><a href=\"#".$interactions."\">Interactions</a></li><li><a href=\"#" . $ora . "\">ORA</a></li><li><a href=\"#$motif\">Motifs</a></li></ul><div id=\"".$network."\"  class=\"innerPanels\">";

  $pb .= "<div class=\"cytoNav\">";

  //$pb .= "<a href=\"#\" onclick=\"$('#".$network."Pane').cytoscape('get').reset();$('#".$network."Pane').cytoscape('get').fit();\">Reset</a>";
  //$pb .= "<a href=\"#\" onclick=\"$('#".$network."Pane').cytoscape('get').center();\">Center</a>";

  if ($dataset == 'disease' || $dataset == 'all') {
    $pb .= "<br><span class=\"padLeft\">Motifs</span>";

      if (!file_exists($subfolder . '/motifs.txt')) {
      $pb .= "<a id=\"startMotifSearchButton\" href=\"#\" title=\"Attention: this may take a while!\" onclick=\"motifSearch('$dataset'); return(false);\">Search for motifs</a><span id=\"motifloader\"><img height=\"15\" src=\"img/loading.gif\"></span>";
    } 
    $pb .= "<select id=\"motifs\" onChange=\"highlightMotif(this.value);\"></select><a href=\"#\" onclick=\"showCos('".session_id()."','".$dataset."','cotargeted'); return(false);\">Show co-targeted</a><a href=\"#\" onclick=\"showCos('".session_id()."','".$dataset."','coregulated'); return(false);\">Show co-regulated</a>";
  }
  
  $pb .= "<br/><span class=\"padLeft\">Layout</span><a href=\"#\" onclick=\"$('#".$network."Pane').cytoscape('get').fit();\">Reset</a>";
  $pb .= "<select name=\"layoutMethod\" onchange=\"$('#".$network."Pane').cytoscape('get').layout({name:$('select[name=layoutMethod] :selected').attr('value')})\">" .getLayouts() . "</select>";
  $pb .= "<a href=\"\" onclick=\"extendSelectedNeighbours(); return false;\">Extend neighbours</a><a href=\"#\" onclick=\"$('#".$network."Pane').cytoscape('get').elements().unselect(); return false; \">Deselect all</a><a href=\"#\" onclick=\"exportImage(); return false;\" title=\"Export current view as PNG image\">Export image</a><br/>";

  if ($dataset =='disease' || $dataset == 'all') {
    $pb .= "<span class=\"padLeft\">Hotspot nodes:</span>";//"<a href=\"#\" onclick=\"return false;\">Hotspot nodes:</a>";
    $pb .= "<a href=\"#\" onclick=\"showDominatingSet('$dataset'); return(false);\">Dominating set</a><span id=\"domsetloader\"><img height=\"15\" src=\"img/loading.gif\"></span>";
  }

  $summary = (parseSummaryFile($subfolder .'/summary.txt')); 

  foreach ($summary as $key => $value) {
    if (stristr($key, 'hotspot') !== FALSE) {

      $hotspots = mergeHotspotsForCytoscape($value);
      $pb .= "<a href='#' onclick=\"$('#".$network."Pane').cytoscape('get').elements().unselect(); $('#".$network."Pane').cytoscape('get').elements('".$hotspots."').select(); return false;\">".getValue($key)."</a>";
      unset($summary[$key]);
    }
  }
  $pb .= "</div>";



  $pb .= "<div id=\"networkParentPane\"><div id=\"".$network."Pane\" class=\"cytoPanel\"></div><div>Selected Nodes:<div id=\"selectedNodes\"></div></div><hr></div><div id=\"summsidebar\">"; 
  
  //$pb .= "<div id=\"".$summary."\"><pre>";

  
  $pb .= "<table>";
  foreach ($summary as $key => $value) {
    //if (preg_match('/(([\w]*(\\|\/){1}){3}[\w]*\.png)/', $value, $matches)) {
    if (preg_match('/(\/.*)*(\/uploads\/([\w-]*\/)*(\w*\.png))/', $value, $matches)) {
      //$pb .= $key . ' ' . $matches[0] . $matches[1] . $matches[2];
      $imgsrc = './uploads/'. session_id() .'/' . $dataset . '/' . $matches[4];
      $pb .= '<tr><td>' . getValue($key) . '</td><td><a class="image-popup-vertical-fit" href="'.$imgsrc.'" title=""><img src="'.$imgsrc.'" height="75"></a></td></tr>';
    } elseif ($value != '') {
     $pb .= '<tr><td>' . getValue($key) . '</td><td>' . urlCheck($value, $dic[$key]) .'</td></tr>';
    }
  }
  $pb .= "</table>";
  
  $pb .= "</div></div>";

  $pb .= "<div id=\"".$interactions."\">";
  $pb .= "<table id=\"".$resTable."\" class=\"stripe\">";
  $pb .= parseTableContentFromCSV($subfolder . '/res.txt');
  $pb .= "</table></div>";

  $pb .= "<div id=\"".$ora."\">";

  if (file_exists($subfolder . '/genes.ora.txt')) {

    $linkTemplates = 
      array(
        'david.functional.clust.link' => 'http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&tool=term2term&annot=GOTERM_BP_ALL&ids=',
        'david.BP.link' => 'http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&tool=chartReport&annot=GOTERM_BP_ALL&ids=',
        'david.KEGG.link' => 'http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&tool=chartReport&annot=KEGG_PATHWAY&ids=',
        'david.OMIM.link' => 'http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&tool=chartReport&annot=OMIM_DISEASE&ids='
        );

    $pb .= "<h4>Gene ORA</h4><h5>External links</h5><table id=\"geneORATable\">";
    $geneOra = (file_get_contents($subfolder .'/genes.ora.txt'));
    $pb .= "<thead><td>Type</td></thead>";
    foreach ($linkTemplates as $key => $value) {
      $pb .= "<tr><td><a target=\"_blank\" href=\"".$value.$geneOra."\">".getValue($key)."</a></td></tr>";
      //$pb .= "<a class=\"externalORA\" href=\"".$value."\" onclick=\"window.open('".$value."', '_blank');\">".getValue($key)."</a><br/>";
    }
    $pb .= "</table>";
  } else {
    $pb .= "Sorry, no data available for gene over representation analysis";
  }

  if (file_exists($subfolder . '/mirna.ora.txt')) {
    $pb .= "<h4>miRNA ORA</h4><table id=\"".$oraTable."\" class=\"stripe\">";
    $pb .= parseTableContentFromCSV($subfolder . '/mirna.ora.txt');
    $pb .= "</table>";
  } else {
    $pb .= "Sorry, no data available for this";
  }

  


  $pb .= "</div>";

  $pb .= "<div id=\"$motif\">";

  if (!file_exists($subfolder . '/motifs.txt')) {
    $pb .= "Motifs have not been computed yet. Please run motif search in network view! Please note this may take a while.";
  } else {
    $pb .= "<table id=\"$motifTable\" class=\"stripe\">";
    $pb .= parseTableContentFromCSV($subfolder."/motifs.txt");
    $pb .= "</table>";
    $pb .= '<script type="text/javascript">motifSearch(\''.$dataset.'\')</script>';

    $pb .= '<script type="text/javascript">';
    $pb .= '$("#'.$dataset.'-motifTable thead").each(function (index) { $(this).find("th").eq(11).after("<th>Co-Regulated</th><th>Co-Targeted</th>")});';
    $pb .= '$("#'.$dataset.'-motifTable tr").each(function(index) {' . "\n";
    $pb .= 'var tf = $(this).find("td").eq(2).text();' . "\n";
    $pb .= 'var mirna = $(this).find("td").eq(3).text();' . "\n";
    $pb .= 'var gene = $(this).find("td").eq(4).text();' . "\n";
    $pb .= '$(this).find("td").eq(11).after("<td><a target=\"_blank\" href=\"motif.php?tf="+tf+"&gene="+gene+"&mirna="+mirna+"&type=coregulated&folder='.session_id().'&dataset='.$dataset.'\">Show co-regulated</a></td><td><a target=\"_blank\" href=\"motif.php?tf="+tf+"&gene="+gene+"&mirna="+mirna+"&type=cotargeted&folder='.session_id().'&dataset='.$dataset.'\">Show co-targeted</a></td>");})';
    $pb .= "</script>";
  }
  $pb .= "</div>";
  $pb .= "</div>";

  //$pb .= "<script type=\"text/javascript\"> $( document ).ready(function() { var cy = $('#".$network."Pane').cytoscape('get');}); </script>";
  return $pb;
}



//<script src="cytoscape.js" type="text/javascript"></script>
//</head>
//<body>
//  <h2>Transcription Factor Prediction</h2>
//  <h3><?php echo getValue($_REQUEST['dataset']); Interactions</h3>
//   
//  <script src="cytoscape.php?folder=<?php echo $_REQUEST['folder'];" type="text/javascript"> </script>

function getLayouts() {

  $layouts = array('concentric','cose','random','circle','grid','breadthfirst');//, 'arbor');

  $pb = '';
  foreach($layouts as $layout) {
    $pb .= "<option value=\"".$layout."\">$layout</option>";
  }

  return $pb;
}

function mergeHotspotsForCytoscape($str) {
  $hotspots = explode(',',$str);

  $pbArr = array();
  foreach ($hotspots as $hotspot) {
    array_push($pbArr, "node[id=\'".trim($hotspot)."\']");
  }

  return implode(',', $pbArr);
}

?>

