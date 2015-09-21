<?php

include 'httpful.phar';

$response = \Httpful\Request::get('http://localhost:1234/v1')->send();

echo $response . "\n";


?>
