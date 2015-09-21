<?php
session_start();
require_once '../Requests/library/Requests.php';
Requests::register_autoloader();

include("createCytoscapeNetwork.php");

$baseDir = realpath("." );
$subfolder = $baseDir . '/uploads/'. session_id() . '/';//. $_REQUEST['dataset'];

//var_dump($_REQUEST);

$filelocation = $subfolder . $_REQUEST['dataset'].'/res.txt';

//echo $filelocation;
$json = createJSON($filelocation);

$json = '{ "data" : { "name" : "TFmiR" } , "elements" : ' . $json . ' };';

//echo "<code>".$json."</code>";

$headers = array('Content-Type' => 'application/json');
$response = Requests::post('http://localhost:1234/v1/networks/motifsearch', $headers, $json, array("timeout" => 1200000));

//echo file_get_contents('res.json');

//$response = Requests::post('http://localhost:1234/v1/networks', $headers, file_get_contents('res.json'));
//echo var_dump($response->body); 

//if ($response->success)
echo $response->body;
//else header("Status: 404 Not Found");
?>