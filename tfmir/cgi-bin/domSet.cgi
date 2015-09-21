#!/usr/bin/perl
print "Content-type: text/html\n\n";
use CGI;

my($q) = new CGI;
#-Dcgi.content_type=$CONTENT_TYPE     -Dcgi.content_length=$CONTENT_LENGTH     -Dcgi.request_method=$REQUEST_METHOD     -Dcgi.query_string=$QUERY_STRING     -Dcgi.server_name=$SERVER_NAME     -Dcgi.server_port=$SERVER_PORT     -Dcgi.script_name=$SCRIPT_NAME     -Dcgi.path_info=$PATH_INFO   

$session = $q->param('session');
$interactions = $q->param('interactions');

#if ($session !~ '/[\w]{32}/') {
#	exit(1);
#}

#if ($interactions !~ '/disease/' || $interactions !~ '/all/') {
#	exit(1);
#}

use Cwd;
my $dir = getcwd;

use Cwd 'abs_path';

$binDir = abs_path(getcwd. '/../../');
$inputFile = abs_path(getcwd . '/../uploads/' . $output . '/' . $session . '/' . $interactions . '/res.txt');
$outputFile = abs_path(getcwd . '/../uploads/' . $output . '/' . $session . '/' . $interactions . '/dominatingSet.txt');
$domlog = abs_path(getcwd . '/../uploads/' . $output . '/' . $session . '/' . $interactions);

#print('#java -jar ' . $binDir . '/DominatingSet.jar ' . $inputFile . ' ' . $outputFile);

# add slash for local
$command = 'java -jar ' . $binDir . '/DominatingSet.jar '  . $inputFile . ' ' . $outputFile. ' 1>'.$domlog.'/program.stdout 2>'.$domlog.'/program.stderr';

my $output = system($command);


print '1';#<html><code>'.$command . $output . '</code></html>';
