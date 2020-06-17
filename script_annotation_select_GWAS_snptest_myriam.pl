#!/usr/bin/perl
use Data::Dumper;
use warnings;
use strict;
use List::MoreUtils qw(uniq);

my $annotation = $ARGV[0] ;
my $line1 ;
my $variant = $ARGV[1];
my $fichier_sortie = $ARGV[3] ;
my $line2 ;
my $line3 ;
my $bp = $ARGV[2] ;
my $p_threshold = $ARGV[4] ; 
my $CHR = $ARGV[5] ; 
 


#création de la table de hash avec toutes les données d'annotation
my %gene_annotation = creer_table_de_hash($annotation) ;


#annotation des snps
my %snp_hash = annotation (\%gene_annotation, $variant, $fichier_sortie);



sub creer_table_de_hash {
	my ($nomdefichier) = @_ ;
	my %gene_annotation ;
	open(FILE1, $nomdefichier) or die "open: $!";
	#print $nomdefichier;
	while (defined($line1 = <FILE1>)) {
		if($. > 1){
			my @line_annotation = split( /\s/, $line1);
			#print "@line_annotation";
			my $chr = $line_annotation[2];
			my $start= $line_annotation[3];
			my $end = $line_annotation[4];
			my $gene_name = $line_annotation[1];
			my $gene_p = $line_annotation[5];
			my $gene_info = $gene_name . " $gene_p";
			$gene_annotation{$chr}{$start}{$end} = $gene_info;
		}
	}
close(FILE1);
#print Dumper(\%gene_annotation);
return %gene_annotation
}



sub annotation {
	my ($reference_gene_annotation, $nomdefichier_2, $out) = @_;
	my %gene_annotation_2 = %{$reference_gene_annotation} ;
	my %snp_hash ;
	my $dif_imp ;
	my $dif_start ;
	my $dif_stop ;
	my $dif_start_o ;
	my $dif_stop_o ;
	my $gene_name ;

	open (FSOR, ">$out");
	print FSOR "rs_id p_GWAS chr pos gene(name p_TWAS dist)\n" ;
	open (FILE2, $nomdefichier_2) or die "open: $!";
	while (defined($line2 = <FILE2>)) {
		if($. > 1){
			$line2 =~ s/\s+/\t/g;
			my @snp = split (/\t/, $line2);
			my $snp_info = $snp[21];
			if ($snp_info >= 0.8 && $snp_info ne "-1"){
				my $snp_chr = $snp[2];
				my $snp_position = $snp[3];
				my $snp_pval = $snp[20];
				if ($snp_pval <= $p_threshold){
					my $snp_name = "$snp[1]" . " $snp_pval";
					chomp $snp_name ;
					#je vais parcourir la table de hash pour comparer les positions début et fin (attention! il faut que le start soit inférieur au end)
					foreach my $clef_start (keys %{$gene_annotation_2{$snp_chr}}){
						my $clef_start_moins = $clef_start - $bp ;
						foreach my $clef_stop (keys %{$gene_annotation_2{$snp_chr}{$clef_start}}) {
							my $clef_stop_plus = $clef_stop + $bp ;
							#je vais compléter la table de hash uniquement si le snp est dans le gène + bp 
							if ($clef_start_moins < $snp_position && $clef_stop_plus > $snp_position) {
								$dif_start_o = ($clef_start - $snp_position);
								$dif_stop_o = ($clef_stop - $snp_position);
								$dif_start = abs($clef_start - $snp_position);
								$dif_stop = abs($clef_stop - $snp_position);
								$gene_name = $gene_annotation_2{$snp_chr}{$clef_start}{$clef_stop} ;							
								#Ensuite je tes l'existence dans la table de hash pour compléter ou créer
								if (exists($snp_hash{$snp_chr}{$snp_position}{$snp_name})) {
									#je teste quelle est la distance la plus courte entre start et stop et ma position
									if ($dif_start_o < 0 && $dif_stop_o > 0) {
										$dif_imp = 0 ;
									}
									else {
										if ($dif_stop < $dif_start){
											$dif_imp = $dif_stop;
										}
										else {
											$dif_imp = $dif_start;
										}
									}
									$snp_hash{$snp_chr}{$snp_position}{$snp_name}{$gene_name} = $dif_imp;
								}
								else {
									if ($dif_start_o < 0 && $dif_stop_o > 0) {
										$dif_imp = 0 ;
									}
									else {
										if ($dif_stop < $dif_start){
											$dif_imp = $dif_stop;
										}
										else {
											$dif_imp = $dif_start;#print "$snp_position\n";
										}
									}
									#print "$gene_annotation_2{$snp_chr}{$clef_start}{$clef_stop}";
									$snp_hash{$snp_chr}{$snp_position}{$snp_name}{$gene_name} = $dif_imp ;
								}
							}
						}
					}
				}
			}
		}
	}
	#je relis le fichier pour que ma sortie soit dans le bon ordre... c'est surement tout pourri :)
	close (FILE2);
	open (FILE3, $nomdefichier_2) or die "open: $!";	
	while (defined($line3 = <FILE3>)) {
		if($. > 1){
			$line3 =~ s/\s+/\t/g;
			my @snp_final = split (/\t/, $line3);
			#print $snp_final[2];
			my @dist_f;
			my @dist_s ;
			my $id_rs = "$snp_final[1]" . " $snp_final[20]";
			my $chr_rs = $snp_final[2] ;
			my $pos_rs = $snp_final[3] ;
			chomp $id_rs;
			#chomp $id_rs ;
			#je teste l'existence d'un gènes pour le SNP position chr
			if (exists ($snp_hash{$chr_rs}{$pos_rs}{$id_rs})) {
				my $gene;
				#je créé un boucle pour récupérer les distances à des gènes pour un SNP
				foreach $gene (keys %{$snp_hash{$chr_rs}{$pos_rs}{$id_rs}}) {
					#je créé ma table de distances
					push @dist_f, $snp_hash{$chr_rs}{$pos_rs}{$id_rs}{$gene};
				}
				#je trie ma table de distances
				@dist_s = sort { $a <=> $b } @dist_f ;
				#je récupère les valeurs unique de distance
				my @dist_u = uniq(@dist_s) ;
				print FSOR "$id_rs $chr_rs $pos_rs ";
				#Une boucle pour sur mes distances pour un SNP et je test si la valeur correspond dans ma table de hash, si oui j'écris "gène et distance" dans le fichier
				for (my $i = 0; $i <= $#dist_u; $i++){			
					my $gene2;
					foreach $gene2 (keys %{$snp_hash{$chr_rs}{$pos_rs}{$id_rs}}) {
						my $dist = $snp_hash{$chr_rs}{$pos_rs}{$id_rs}{$gene2};
						if ($dist == $dist_u[$i]){
							print FSOR  "$gene2  $dist " ;
						}
					}
		
				}
		
			print FSOR "\n";
			}
		}
	}
close (FILE3);
close (FSOR);
#print Dumper(\%snp_hash);
return %snp_hash ;
}

