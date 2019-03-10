use strict;
package fye;

# Routines for computing the Fisher-Yates exact test in the log-domain.

use constant PI => 4*atan2(1,1);

# lf(n) := ln(fact(n)) , according to Ramanujan's approximation.
sub lf {
	my $n = shift(@_);
	return ($n > 1) ? $n*log($n) - $n + log($n*(1 + 4*$n*(1 + 2*$n)))/6 + log(PI)/2 : 0;
}

# Definition of Fisher-Yates exact test for 2x2 contingency table ((a b) (c d)).
sub lnFYabcd {	
	my ($a, $b, $c, $d) = @_;
	return lf($a+$b)-lf($a) + lf($c+$d)-lf($d) + lf($a+$c)-lf($c) + lf($b+$d)-lf($b) - lf($a+$b+$c+$d);
}

# With Evert's notation: O11 = a; R1 = a+b; C1 = a+c; N = a+b+c+d  ; a = O11; b = R1-O11; c = C1-O11; d=N-R1-C1+O11
sub lnFY {	
	my ($O11, $R1, $C1, $N) = @_;
	return lf($R1)-lf($O11) + lf($N-$R1)-lf($N-$R1-$C1+$O11) + lf($C1)-lf($C1-$O11) + lf($N-$C1)-lf($R1-$O11) - lf($N);
}

# Calculate (approximately) the tail log-probabilities.
sub lnFYtail {
	my ($O11, $R1, $C1, $N) = @_;
	my $lnP = lnFY($O11, $R1, $C1, $N);
	my $dir = ($O11 < $R1*$C1/$N) ? 1 : -1;

	# Update the tail with the log-sum-exp trick (should use log1p() if available).
	if ($dir > 0) {
		$lnP += log(1 + exp(lnFY($O11-1, $R1, $C1, $N) - $lnP)) if $O11 >= 1;
		$lnP += log(1 + exp(lnFY($O11-2, $R1, $C1, $N) - $lnP)) if $O11 >= 2;
		$lnP += log(1 + exp(lnFY($O11-3, $R1, $C1, $N) - $lnP)) if $O11 >= 3;
		$lnP += log(1 + exp(lnFY($O11-4, $R1, $C1, $N) - $lnP)) if $O11 >= 4;
		$lnP += log(1 + exp(lnFY($O11-5, $R1, $C1, $N) - $lnP)) if $O11 >= 5;
	} else {
		$lnP += log(1 + exp(lnFY($O11+1, $R1, $C1, $N) - $lnP)) if $C1 - $O11 >= 1;
		$lnP += log(1 + exp(lnFY($O11+2, $R1, $C1, $N) - $lnP)) if $C1 - $O11 >= 2;
		$lnP += log(1 + exp(lnFY($O11+3, $R1, $C1, $N) - $lnP)) if $C1 - $O11 >= 3;
		$lnP += log(1 + exp(lnFY($O11+4, $R1, $C1, $N) - $lnP)) if $C1 - $O11 >= 4;
		$lnP += log(1 + exp(lnFY($O11+5, $R1, $C1, $N) - $lnP)) if $C1 - $O11 >= 5;
	}
	return ($lnP, $dir);
}

# Minimum Sensitivity
sub minSens {	
	my ($O11, $R1, $C1, $N) = @_;
	return ($R1 > $C1) ? $O11/$R1 : $O11/$C1;
}

# Package success!
1;
