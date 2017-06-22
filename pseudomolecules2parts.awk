#!/usr/bin/mawk -f

BEGIN {
  FS=OFS="\t"
}
{
  if(NR==FNR){ #FIRST INPUT FILE
    header=header"##sequence-region "$1" 1 "$2"\n";
    if($1 ~ /_part1/){sub("_part1","");
      len[$1]=$2
    }
  } else { #SECOND INPUT FILE
    if(FNR==2) {
      printf header
    }
    if($1 ~ /^##sequence-region/ ) 
    {
      ##DON'T PRINT
    } else if($1 ~/^#/ ) 
    {
      print ##OTHER COMMENT LINES
    } else if($1 in len && $4>len[$1]) { 
      $4-=len[$1]+1;
      $5-=len[$1]+1;
      $1=$1"_part2";
      print
    } else {
      if($1 in len) {
        $1=$1"_part1"
      }
      print
    }
  }
}
