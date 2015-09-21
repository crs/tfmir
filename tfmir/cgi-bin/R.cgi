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

#%cookies = CGI::Cookie->fetch;
#$miRNA = $cookies{'miRNA'}->value;
#$tf = $cookies{'mRNA'}->value;

# Reconstruct any inputs
#if ($results=$q->param('CALLER')) {
#    $results='<INPUT TYPE=HIDDEN NAME="CALLER" VALUE="'.$results.'">';
#    }
#@input = $q->param('INPUT');

#$in=join("\n",$q->param('INPUT'));

#$in=$q->param('INPUT');
#TFMir(tf.path,mirna.path,pval.cutoff,disease,output.folder)

#$miRNA = $q->param('miRNA-filename');

#$pval = $q->param('orapvalue');
#$disease = $q->param('disease');
#$output = $q->param('session');


###### FIX THIS!!! #####
$temp = join('', "/Users/chris/Website/TF-site/tf/tmp", $output);
#mkdir($temp, 0777);

#$output = '../tf/uploads/' . $output;
#$arguments = join("','", $tf, $miRNA, $pval, $disease, $output);

$in="getwd();";#"source('/Users/chris/Website/TF-site/tf/src/interface.R'); start('".$arguments."');";

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

my $x = system(`R --vanilla -q <$temp/R3.Rin >$temp/R3.Rout`);

# Read the output
    open(IN, "$temp/R3.Rout");

# Tidy it up for HTML
    $results =  join(" ",<IN>);
    $results =~ s/&/&amp;/g;
    $results =~ s/</&lt;/g;
    $results =~ s/>/&gt;/g;
    close(IN);

}
