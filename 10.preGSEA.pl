use strict;
#use warnings;

######Usage: perl preGSEA.pl VCAN

my $colNum=0;
my $rowNum=0;
my %hash=();
my $geneName=$ARGV[0];
my @indexs=();
my @geneArr=();

open(RF,"symbol.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	my @arr=split(/\t/,$line);
	if($.==1){
		for(my $i=1;$i<=$#arr;$i++){
			my @samples=split(/\-/,$arr[$i]);
			if($samples[3]=~/^0/){
			  push(@indexs,$i);
			  my $sampleName=$arr[$i];
			  $hash{$sampleName}=1;
			  $colNum++;
			}
		}
	}
	else{
			$rowNum++;
		  if($arr[0] eq $geneName){
			  foreach my $col(@indexs){
				  push(@geneArr,$arr[$col]);
			  }
		  }
	}
}
close(RF);

my $firstGeneVal=$geneArr[0];
my $geneMed=median(@geneArr);
open(RF,"symbol.txt") or die $!;
open(WF,">$geneName.gct") or die $!;
print WF "#1.2\n";
print WF "$rowNum\t$colNum\n";
open(CLS,">$geneName.cls") or die $!;
print CLS "$colNum\t2\t1\n";

if($firstGeneVal>$geneMed){
	print CLS "#\th\tl\n";
}
else{
	print CLS "#\tl\th\n";
}

@indexs=();
my @typeArr=();
while(my $line=<RF>){
	chomp($line);
		my @arr=split(/\t/,$line);
	if($.==1){
		print WF "NAME\tDESCRIPTION";
		for(my $i=1;$i<=$#arr;$i++){
			my @samples=split(/\-/,$arr[$i]);
			if($samples[3]=~/^0/){
			  my $sampleName=$arr[$i];
			  if(exists $hash{$sampleName}){
				  push(@indexs,$i);
				  print WF "\t$arr[$i]";
				  #delete($hash{$sampleName});
			  }
			}
		}
		print WF "\n";
	}
	else{
		my $symbolName=$arr[0];
		$symbolName=~s/(.+?)\|.+/$1/g;
		print WF "$symbolName\tna";
		foreach my $col(@indexs){
			print WF "\t$arr[$col]";
		}
			print WF "\n";
		if($arr[0] eq $geneName){
			foreach my $col(@indexs){
				if($arr[$col]>$geneMed){
				  push(@typeArr,"h");
			  }
			  else{
			  	push(@typeArr,"l");
			  }
			}
		}
	}
}
print CLS join("\t",@typeArr) . "\n";
close(WF);
close(CLS);
close(RF);

sub median{  
    my (@data) = sort {$a <=> $b} @_;  #  ԭ              
    if(scalar (@data) % 2){  
        return ($data [@data / 2]);  #@data/2 ķ   ֵ    ȡ    
    }else {  
        my ($upper ,$lower);  
        $upper = $data[@data / 2];  
        $lower = $data[@data / 2 -1];  
        return (($lower+$upper) / 2);  
    }  
}  


######      ѧ  : https://www.biowolf.cn/
###### γ     1: https://shop119322454.taobao.com
###### γ     2: https://ke.biowolf.cn
###### γ     3: https://ke.biowolf.cn/mobile
###### ⿡  ʦ   䣺seqbio@foxmail.com
###### ⿡  ʦ΢  : seqBio
