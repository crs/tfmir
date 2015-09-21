#!/usr/bin/perl

# hello.pl -- my first perl script!

print "Content-type: text/html\n\n";

print "Hello, world!\n";

#my $x = `R  --vanilla -e 'print(version)'`;
#my $x = `ls`;
#my $x = `ps -ef | grep dropbox`;
my $x = `/usr/bin/R -e "print(version)" > stupid.txt`;
#my $x = `. /Users/chris/.profile && /usr/local/bin/R -e 'print(version)'`;

print($x)
