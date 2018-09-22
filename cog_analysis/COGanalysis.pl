#!usr/bin/perl
use 5.010;

# Prepare a list of file names
@files = qw(
Ensifer_adhaerens_Casida_A.faa
Ensifer_adhaerens_OV14.faa
Sinorhizobium_sp_A49.faa
Sinorhizobium_sp_M14.faa
);

# Array of each COG category
@COGs = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z);

# Run the analysis
foreach $i (@files) {
	$total = 0;
	$input = "Proteomes/$i";
	open($in, '<', $input);
	while(<$in>) {
		if(/>/) {
			$total++;
		}
	}
	close($in);
	$i =~ s/.faa//g;
	system("cut -f11 emapperOutput/$i.emapper.annotations | tail -n +4 | head -n -3 > COG_categories/$i.txt");
	foreach $j (@COGs) {
		$COGcount = 0;
		$input = "COG_categories/$i";
		$input = $input . '.txt';
		open($in, '<', $input);
		while(<$in>) {
			if(/$j/) {
				$COGcount++;
			}
		}
		close($in);
		$output = "COG_categories/$i.counted.txt";
		$percent = 100 * $COGcount / $total;
		open($out, '>>', $output);
		print $out ("$j\t$COGcount\t$total\t$percent\n");
		close($out);
	}
}

