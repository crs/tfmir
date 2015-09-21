<?php
//echo "<pre>";
//echo "</pre>";
function createJSON($filename) {

	$nodeTemplate = "\n" .'{ "data": { "id": "%s", "type":"%s", "in_disease":"%s", "regulation":"%s", "shape":"%s", "color":"%s"} }';
	$edgeTemplate = "\n" .'{ "data": { "id": "%d", "source":"%s", "target":"%s", "interaction":"%s", "evidence":"%s", "dataSource":"%s", "name":"%s (%s) %s" } }';

  //$nodeTemplate = "\n{ \"data\": { \"id\": \"%s\"', type:'%s', in_disease:'%s', regulation:'%s', shape:'%s', color:'%s'} }";
  //$edgeTemplate = "\n{ \"data\": { \"id\": '%d', source:'%s', target:'%s', type:'%s', evidence:'%s', dataSource:'%s' } }";


	if (!file_exists($filename) || ($f = fopen($filename, "r")) === false) { return "Not available for this type of interaction"; };
  
  	// get table header
  	$first = fgetcsv($f, 0, "\t");

  	$header = array();
  	foreach($first as $head) array_push($header, $head);

	$nodes = array();

	$edges = array();
	$edge_id = 0;
	while (($line = fgetcsv($f, 0, "\t")) !== false) {
  		$cell_no = 0;
  		//foreach($line as $cell) {
  		//	if ($header[$cell_no] == 'regulator' | $header[$cell_no] == 'target') {
  		//		array_push($nodes, $cell)
  		//	}
  		
  		$node_1 = $line[array_search('regulator', $header)];
  		$node_2 = $line[array_search('target', $header)];

  		$categoryString = $line[array_search('category', $header)];

  		$categories = explode('-',$categoryString);

  		$reg_in_disease = $line[array_search('is_regulator_in_disease', $header)];
  		$target_in_disease = $line[array_search('is_target_in_disease', $header)];

  		if (array_search('regulator.reg', $header) !== false) {
  			$regulator_regulation = $line[array_search('regulator.reg', $header)];
  			$target_regulation = $line[array_search('target.reg', $header)];
  		} else {
  			$regulator_regulation = 0;
  			$target_regulation = 0;
  		}

  		$evidence = $line[array_search('evidence', $header)];
  		$dataSource = $line[array_search('source', $header)];

  		$categoryShape = array('tf' => 'hexagon', 'mirna' => 'ellipse', 'gene' => 'rectangle' );
  		$categoryColor = array('tf' => '#66d15c', 'mirna' => '#fffb7d', 'gene' => '#f2780c' );

  		array_push($nodes, sprintf($nodeTemplate, $node_1, $categories[0], $reg_in_disease, $regulator_regulation, $categoryShape[$categories[0]], $categoryColor[$categories[0]]));
  		array_push($nodes, sprintf($nodeTemplate, $node_2, $categories[1], $target_in_disease, $target_regulation, $categoryShape[$categories[1]], $categoryColor[$categories[1]]));  		
  		array_push($edges, sprintf($edgeTemplate, $edge_id, $node_1, $node_2, $categoryString, $evidence, $dataSource, $node_1, $categoryString, $node_2));

//  		$categories = explode('-',$line[2]);
//  		array_push($nodes, sprintf($nodeTemplate, $line[0], $categories[0], $line[7], $line[9]));
//  		array_push($nodes, sprintf($nodeTemplate, $line[1], $categories[1], $line[8], $line[10]));  		
//  		array_push($edges, sprintf($edgeTemplate, $edge_id, $line[0], $line[1], $line[2], $line[3], $line[4]));
  		$edge_id += 1;
  	}
  	$nodes = array_unique($nodes);
  	//echo '{ "data": { "name":"TFmir" },'. "\n";
    //echo '"elements": {';

    $pb =  "{";
    $pb .= '"nodes" : [ ' . implode($nodes, ",") . "],\n";
	  $pb .= '"edges" : [ ' . implode($edges, ",") . "]";
    $pb .= "}";
    return($pb);
    //echo "}\n };";
}

?>