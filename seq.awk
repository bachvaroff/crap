#!/usr/bin/awk -f

BEGIN {
	step = 1;
	
	if (ARGC < 2) {
		print "no args"
		exit;
	} else if (ARGC == 2) {
		from = 1;
		to = ARGV[1];
	} else if (ARGC == 3) {
		from = ARGV[1];
		to = ARGV[2];
	} else if (ARGC >= 4) {
		from = ARGV[1];
		step = ARGV[2];
		to = ARGV[3];
	}
	
	if (!step) step = 1;

	for (i = from; i <= to; i += step)
		print i;
}
