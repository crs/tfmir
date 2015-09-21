#!/usr/bin/perl
#
# Rcgi main file...  see documentation for licencing.

my($temp);

#$temp = "/mnt/disknt1/statguest/statweb/public_html/tmp";
#$temp = "/opt/lampp/htdocs/tf/tmp";


$allowedhost = '';
# $allowedhost = '\.edu';


#use strict;    # the referer() seems to cause strict problems
use CGI;
use CGI::Cookie;
my($results, $in);
my($q)= new CGI;

%cookies = CGI::Cookie->fetch;

my($miRNA);
my($mRNA);

$miRNA = $cookies{'miRNA'}->value;
$mRNA = $cookies{'mRNA'}->value;

# Reconstruct any inputs
#if ($results=$q->param('CALLER')) {
#    $results='<INPUT TYPE=HIDDEN NAME="CALLER" VALUE="'.$results.'">';
#    }
#@input = $q->param('INPUT');

#$in=join("\n",$q->param('INPUT'));

#$in=$q->param('INPUT');
#TFMir(tf.path,mirna.path,pval.cutoff,disease,output.folder)

#$miRNA = $q->param('miRNA-filename');

$pval = $q->param('orapvalue');
$disease = $q->param('disease');
$output = $q->param('session');
$evidence = $q->param('evidence');


###### FIX THIS!!! #####


use Cwd;
my $dir = getcwd;

use Cwd 'abs_path';

$temp = abs_path(getcwd . '/../uploads/' . $output);
mkdir($temp, 0777);

$mRNA = abs_path(getcwd . '/../uploads/' . $output . '/mRNA/' . $mRNA);
$miRNA = abs_path(getcwd . '/../uploads/' . $output . '/miRNA/' . $miRNA);

$output = abs_path( getcwd . '/../uploads/' . $output );
$arguments = join("','", $mRNA, $miRNA, $pval, $evidence, $disease, $output);

my $importPath = abs_path(getcwd . "/../src/interface.R");
$in="source('" . $importPath . "'); start('".$arguments."');";


# Print input form
print $q->header,
'<html>';
#<B>Program Input</B><BR>
#<FORM METHOD=POST ACTION="'.$q->url.'">
#<TEXTAREA NAME="INPUT" ROWS=5 COLS=64>
#'.$in.'</TEXTAREA><BR><INPUT TYPE=SUBMIT VALUE="Go!">
#'.$results.'</FORM>';

# Basic security checks
if (($q->referer()) && ($q->referer() !~ /$allowedhost\//)) {
  print '<b>Sorry, the website which sent you here is not allowed
to run scripts here.  Email the author of the previous page and ask them
to look into it.</b>';
exit;
}

if ($in) {
# Give our input (if any) to R
    chdir($temp);
#    open(R, "| tee $temp/R3.Rin | R --no-save >$temp/R3.Rout");
    open(R, ">$temp/R3.Rin");
    $in =~ s/\r//g;
    $in =~ s/(unlink|postscript|bitmap|unix|system|pdf|pictex|xfig)[^\n]*//g;
    $in =~ s/(file|dir.create|cat|open|close|pipe|fifo|socket)[^\n]*//g;
    # I think they're the only really nasty things to trap...
    print R $in."\n"; # End it politely (not an EOF at the end of a command)
    close R;

my $x = system(`/usr/local/bin/R --vanilla -q <$temp/R3.Rin >$temp/R3.Rout`);

# Read the output
    open(IN, "$temp/R3.Rout");

# Tidy it up for HTML
    $results =  join(" ",<IN>);
    $results =~ s/&/&amp;/g;
    $results =~ s/</&lt;/g;
    $results =~ s/>/&gt;/g;
    close(IN);

# Print it out
    print $temp .'<br>'. $output . '\n' . $mRNA . '\n' . $miRNA . '\nSubmitted job '.$$.':'.`date +%Y%m%d%H%M%S`;#.$results;

# Have we got a call-back name?
#    if ($in=$q->param('CALLER')) {
#        print '<HR><A HREF="'.$in.'">Return to '.$in.'</A>';
#    }
   }

# print '<HR><A HREF="/Rdoc/doc/html/">R language help</A></BODY></HTML>';
 print '</html>';
