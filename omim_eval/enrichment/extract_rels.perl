#!/usr/bin/perl
use strict; use warnings; use autodie;

open my $handle, '<', qq(hp.obo);
chomp(my @lines = <$handle>);
close $handle;

foreach (@lines) {

	    if ($_=~ /^id/) {
	        print(qq($_\n));
	    }
	        if ($_=~ /^name/) {
	        print(qq($_\n));
	    }
	        if ($_=~ /^is_a/) {
	        $_ =~ s/ !.*//g;
	        print(qq($_\n));
	    }
	}
